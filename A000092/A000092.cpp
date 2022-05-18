#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdio>
#include <iostream>
#include <numeric>
#include <vector>

using std::vector;

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

    // Can build sorted list of (j, k) pairs to better help with locality
    uint32_t max_i = isqrt(N);
    vector<uint32_t> pairs;
    for (uint32_t i = 1; i <= max_i; i++) {
        if (i % 100 == 0) {
            std::cerr << i << "/" << max_i << "  " << pairs.size() << std::endl;
        }
        uint32_t i_2 = i*i;

        // Remove anything where i_2 + pair > N
        while (!pairs.empty() && pairs.back() + i_2 > N) {
            pairs.pop_back();
        }

        // add new pairs for j = i-1
        {
            uint32_t j = i - 1;
            uint32_t j_2 = j*j;
            uint32_t i_j = i_2 + j_2;
            if (j > 1 && i_j + 1 <= N) {
              size_t start = pairs.size();

              for (uint32_t k = 1; k < j; k++) {
                  uint32_t k_2 = k*k;
                  if (i_j + k_2 > N) break;
                  pairs.push_back(j_2 + k_2);
              }
              assert(pairs.size() > start);

              auto merge_start = std::lower_bound(pairs.begin(), pairs.begin() + start, j_2 + 1);
              assert(abs(std::distance(pairs.begin(), merge_start)) <= start);
              std::inplace_merge(merge_start, pairs.begin() + start, pairs.end());
              //assert(std::is_sorted(pairs.begin(), pairs.end()));
            }
        }

        tuples += pairs.size();

        // Data being sorted means access to counts is sequential and fast
        for (auto p : pairs) {
            uint32_t n = i_2 + p;
            assert(n <= N);
            counts[n] += 48;  // 6 * 8
        }
    }

    printf("\ttuples: %lu\n", tuples);
    //printf("\tsum(counts): %lu\n", sum(counts))
    //printf("\tmax(counts): %lu\n", max(counts))

    return counts;
}

void enumerate_n3(uint32_t N) {
    auto counts = get_n3_counts_v2(N);

    if (N > 1000) {
        // Nice verification check from A117609
        uint32_t t = 0;
        for (uint32_t i = 0; i <= 1000; i++) t += counts[i];
        std::cerr << "\tSum of 0..1000: " << t << std::endl;
        assert(t == 132451);
    }

    printf("| nth | n = A000092 | P(n) = A000223 | A(n) = A000413 |\n");

    uint64_t A_n = 0;
    double record = 1;  // To avoid initial zero term
    vector<uint32_t> A000092;
    vector<uint32_t> A000223;
    vector<uint32_t> A000413;
    for (uint32_t n = 0; n <= N; n++) {
        A_n += counts[n];

        double V_n = 4.0/3.0 * M_PI * pow(n, 1.5);
        double P_n = fabs(A_n - V_n);
        if (P_n > record) {
            A000092.push_back(n);
            double P_n_rounded = round(P_n);
            A000223.push_back(P_n_rounded);
            A000413.push_back(A_n);
            record = P_n;
            uint32_t nth = A000092.size();
            if ((nth < 10) || (nth % 5 == 0) || (nth > 120)) {
                printf("| %3d | %11d | %14.0f | %14lu |\n", nth, n, P_n_rounded, A_n);
            }
        }
    }

//    for fn, An in [("b000092.txt", A000092), ("b000223.txt", A000223), ("b000413.txt", A000413)]:
//        with open(fn, "w") as f:
//            for i, a in enumerate(An, 1):
//                f.write(f"{i} {a}\n")
//
}

int main(void) {
    // For 100 terms in 1 second
    //enumerate_n3(1560000);
    // For 124 terms in 34 seconds
    //enumerate_n3(10000000);

    enumerate_n3(45 * 10'000'000);
}
