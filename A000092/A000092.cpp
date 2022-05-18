#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdio>
#include <iostream>
#include <numeric>
#include <vector>
#include <set>

using std::vector;
using std::multiset;

uint32_t isqrt(uint32_t n) {
    uint32_t t = sqrt(n);
    while (t*t > n) {
      t -= 1;
    }
    while ((t+1) * (t+1) <= n) {
      t += 1;
    }
    return t;
}

uint32_t* get_n3_counts_v2(size_t N) {
    // Could store counts as uint16_t counts_48 + uint16_t counts * 4
    // would make counts_48 access twice as local
    uint32_t* counts = (uint32_t*) malloc(sizeof(uint32_t) * (N+1));
    std::fill(counts, counts + N+1, 0);

    size_t tuples = 0;
    // i == j == k
    counts[0] += 1;
    for (uint32_t i = 1; i <= isqrt(N / 3); i++) {
        uint32_t n = 3*i*i;
        assert(n <= N);
        tuples += 1;
        counts[n] += 8;
    }

    // i = j, j > k
    for (uint32_t i = 1; i <= sqrt(N / 2); i++) {
        uint32_t temp = 2*i*i;
        assert(temp <= N);

        // k = 0
        tuples += 1;
        counts[temp] += 3 * 4;

        // k > 0, k < j
        for (uint32_t k = 1; k < i; k++) {
            uint32_t n = temp + k*k;
            if (n > N) break;
            tuples += 1;
            counts[n] += 24;  // 3 * 8
        }
    }

    // i > j = k
    for (uint32_t i = 1; i <= isqrt(N); i++) {
        uint32_t temp = i*i;
        assert(temp <= N);

        // j = k = 0
        tuples += 1;
        counts[temp] += 6;  // 3 * 2

        // j = k, j > 0
        for (uint32_t j = 1; j < i; j++) {
            uint32_t n = temp + 2*j*j;
            if (n > N) break;
            tuples += 1;
            counts[n] += 24;  // 3 * 8
        }
    }

    for (uint32_t i = 1; i <= isqrt(N); i++) {
        uint32_t i_2 = i*i;
        // i > j, k = 0
        for (uint32_t j = 1; j < i; j++) {
            uint32_t i_j = i_2 + j*j;
            if (i_j > N)
                break;
            // k = 0
            tuples += 1;
            counts[i_j] += 24;  // 6 * 4
        }
    }

    uint32_t max_i = isqrt(N);

    /**
     * Build sorted list of (j_2 + k_2) to better help with locality
     * `fixed_pairs` contains all entries <= j*j
     * `merging` contains elements larger
     */
    vector<uint32_t> fixed_pairs;
    multiset<uint32_t> merging;

    for (uint32_t i = 1; i <= max_i; i++) {
        uint32_t i_2 = i*i;

        size_t popped = 0, fixed_pos = 0, pushed = 0, merged = 0;
        // Remove anything where i_2 + pair > N
        while (!fixed_pairs.empty() && fixed_pairs.back() + i_2 > N) {
            fixed_pairs.pop_back();
            popped++;
        }
        assert(merging.empty() || *merging.rbegin() <= N);
        merging.erase(merging.upper_bound(N - i_2), merging.end());
        assert(merging.empty() || (*merging.rbegin() + i_2) <= N);

        // add new pairs for j = i-1
        {
            uint32_t j = i - 1;
            uint32_t j_2 = j*j;
            uint32_t i_j = i_2 + j_2;
            if (j > 1 && i_j + 1 <= N) {  // TODO test if j > 1 is still needed with asserts
                /*
                size_t start = pairs.size();

                for (uint32_t k = 1; k < j; k++) {
                    uint32_t k_2 = k*k;
                    if (i_j + k_2 > N) break;
                    pairs.push_back(j_2 + k_2);
                }
                assert(pairs.size() > start);
                pushed = pairs.size() - start;

                // Early part of list is constant forever
                auto merge_start = std::lower_bound(pairs.begin(), pairs.begin() + start, j_2 + 1);
                fixed_pos = std::distance(pairs.begin(), merge_start);
                merged = start - fixed_pos;
                assert(abs(std::distance(pairs.begin(), merge_start)) <= start);

                // Merging a few new entries (X00 - X000) into LONG list (> million)
                // Nice that it's O(n) but can it be faster -> binary tree?
                std::inplace_merge(merge_start, pairs.begin() + start, pairs.end());
                */

                //assert(std::is_sorted(pairs.begin(), pairs.end()));
                // Move small items out of merging
                {
                    auto it = merging.begin();
                    for (; it != merging.end() && *it <= j_2; it++) {
                        fixed_pairs.push_back(*it);
                    }
                    merging.erase(merging.begin(), it);
                }

                fixed_pos = fixed_pairs.size();
                merged = merging.size();

                for (uint32_t k = 1; k < j; k++) {
                    uint32_t k_2 = k*k;
                    if (i_j + k_2 > N) break;
                    merging.insert(j_2 + k_2);
                    pushed++;
                }
            }
        }

        if (i && (20 * i) % max_i < 20) {
            fprintf(stderr, "\t%6d/%d pairs: %lu (%lu removed, %lu fixed, %lu new, %lu merged)\n",
                    i, max_i, fixed_pairs.size(), popped, fixed_pos, pushed, merged);
        }

        tuples += fixed_pairs.size() + merging.size();

        // Data being sorted means access to counts is sequential and fast
        for (auto p : fixed_pairs) {
            uint32_t n = i_2 + p;
            assert(n <= N);
            counts[n] += 48;  // 6 * 8
        }
        for (auto p : merging) {
            uint32_t n = i_2 + p;
            assert(n <= N);
            counts[n] += 48;  // 6 * 8
        }
    }

    printf("\ttuples: %lu\n", tuples);

    return counts;
}

void enumerate_n3(uint64_t N) {
    assert(N <= std::numeric_limits<uint32_t>::max());

    auto counts = get_n3_counts_v2(N);

    if (N > 1000) {
        // Nice verification check from A117609
        uint64_t t = 0;
        for (uint32_t i = 0; i <= 1000; i++) t += counts[i];
        std::cerr << "\tSum of 0..1000: " << t << std::endl;
        assert(t == 132451);

        for (uint32_t i = 1001; i <= N; i++) t += counts[i];
        std::cerr << "\tSum(counts)   : " << t << std::endl;
    }

    printf("| nth | n = A000092 | P(n) = A000223 | A(n) = A000413  |\n");

    uint64_t A_n = 0;
    double record = 2;  // To avoid initial zero term
    vector<uint32_t> A000092;
    vector<uint32_t> A000223;
    vector<uint32_t> A000413;
    for (uint32_t n = 0; n <= N; n++) {
        A_n += counts[n];

        // This starts to have rounding error at some point
        double V_n = 4.0/3.0 * M_PI * pow(n, 1.5);
        double P_n = fabs(A_n - V_n);
        if (P_n > record) {
            A000092.push_back(n);
            // Can easily round incorrectly after 200 million
            double P_n_rounded = round(P_n);
            A000223.push_back(P_n_rounded);
            A000413.push_back(A_n);
            record = P_n;
            uint32_t nth = A000092.size();
            if ((nth < 10) || (nth % 5 == 0) || (nth > 120)) {
                printf("| %3d | %11d | %14.0f | %15lu | -> %.5f\n", nth, n, P_n_rounded, A_n, record);
            }
        } else if (P_n + 0.01 > record) {
            printf("near_miss for record | --- | %11d | --- | %15lu | -> %.5f\n", n, A_n, record);

        }
    }

//    for fn, An in [("b000092.txt", A000092), ("b000223.txt", A000223), ("b000413.txt", A000413)]:
//        with open(fn, "w") as f:
//            for i, a in enumerate(An, 1):
//                f.write(f"{i} {a}\n")
//
}

int main(void) {
    uint64_t ONE_MILLION = 1'000'000;

    // For 100 terms in 0.14 second
    //enumerate_n3(1560000);

    // For 124 terms in 3.2 seconds
    enumerate_n3(10 * ONE_MILLION);

    // For 151 terms in 76 seconds
    //enumerate_n3(63 * ONE_MILLION);

    // For 188 terms in 21 minutes
    //enumerate_n3(450 * ONE_MILLION);
}
