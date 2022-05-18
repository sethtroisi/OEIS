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
    assert(t*t <= n);
    assert(t*t + 2*t + 1 > n);
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

    // Improve memory access by using smaller counts
    uint16_t* counts_temp = (uint16_t*) malloc(sizeof(uint16_t) * (N+1));
    std::fill(counts_temp, counts_temp + N+1, 0);
    uint64_t min_temp = 0;

    uint32_t max_i = isqrt(N);
    // Build sorted list of (j^2 + k^2) to better help with locality
    vector<uint32_t> pairs;
    for (uint32_t i = 3; i <= max_i; i++) {
        uint32_t i_2 = i*i;
        uint32_t max_pair = N - i_2;

        size_t popped = pairs.size(), fixed_pos = 0, pushed = 0, merged = 0;
        while (!pairs.empty() && pairs.back() > max_pair) {
            pairs.pop_back();
        }
        popped -= pairs.size();

        // add new pairs for j = i-1
        {
            uint32_t j = i - 1;
            uint32_t j_2 = j*j;
            if (i_2 + j_2 + 1 <= N) {
              size_t start_pairs = pairs.size();

              for (uint32_t k = 1; k < j; k++) {
                  uint32_t pair = j_2 + k*k;
                  if (pair > max_pair) break;
                  pairs.push_back(pair);
              }
              assert(pairs.size() > start_pairs);
              pushed = pairs.size() - start_pairs;
              assert(pushed > 0);

              // Early part of list is constant forever
              {
                  //auto merge_start = std::lower_bound(pairs.begin(), pairs.begin() + start_pairs, j_2 + 1);
                  auto merge_start = std::upper_bound(pairs.begin(), pairs.begin() + start_pairs, j_2 + 1);
                  fixed_pos = std::distance(pairs.begin(), merge_start);
                  assert(fixed_pos <= start_pairs);
                  merged = start_pairs - fixed_pos;
              }

              /**
               * This is merging a few new entries (X00 - X000) into a very LONG list (> million)
               * Nice that it's O(n) but can it be faster?
               * I tried with binary tree hoping that O(small * log(n)) < O(n); way slower out.
               * I tired storing up storing up temp_pairs before merging into pairs; same speed.
               */
              std::inplace_merge(pairs.begin() + fixed_pos, pairs.begin() + start_pairs, pairs.end());
              //assert(std::is_sorted(pairs.begin(), pairs.end()));
            }
        }

        if (i && (20 * i) % max_i < 20) {
            fprintf(stderr, "\t%6d/%d pairs: %lu (%lu removed, %lu fixed, %lu new, %lu merged, %lu processed)\n",
                    i, max_i, pairs.size(), popped, fixed_pos, pushed, merged, tuples);
        }

        tuples += pairs.size();

        // Data being sorted means access to counts is sequential and fast
        // Can almost be parellel here but would need to account for duplicates in pairs
        for (auto p : pairs) {
            uint32_t n = i_2 + p;
            assert(n <= N);
            //counts[n] += 48;  // 6 * 8
            counts_temp[n] += 1;
        }

        /**
         * Numbers may have many representations r3 (https://oeis.org/A025436)
         * Make sure to clear out temp before counts_temp > 65535/4 (or I'll worry about overflow)
         * TODO: could check that sum(counts_temp) = delta tuples
         */
        if (i % 2048 == 0 || i == max_i) {
            uint64_t max_temp = std::min<uint64_t>(3ul * i_2, N+1);
            uint16_t max_seen = *std::max_element(counts_temp + min_temp, counts_temp + max_temp);
            for (uint32_t i = min_temp; i < max_temp ; i++) {
                counts[i] += (uint32_t) 48 * counts_temp[i];
            }
            std::fill(counts_temp + min_temp, counts_temp + max_temp, 0);
            printf("\t\tAt i=%d, cleared %lu, max_seen=%d\n", i, max_temp - min_temp, max_seen);
            min_temp = i_2; // technically (i+1) * (i+1)
            assert(max_seen < 0x3FFF); // Make sure we don't overflow
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
            // Can round incorrectly for large n, use verification.py to validate A0002223
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

    // For 124 terms in 2 seconds
    //enumerate_n3(10 * ONE_MILLION);

    // For 151 terms in 42 seconds
    //enumerate_n3(63 * ONE_MILLION);

    // For 175 terms in ~4 minute test
    enumerate_n3(264 * ONE_MILLION);

    // For 188 terms in 13 minutes
    //enumerate_n3(450 * ONE_MILLION);

    // For 200 terms in XX minutes
    //enumerate_n3(860 * ONE_MILLION);

    // For 210 terms in XX minutes
    //enumerate_n3(1400 * ONE_MILLION);

    // For XXX terms in XX minutes
    //enumerate_n3(2000 * ONE_MILLION);
}
