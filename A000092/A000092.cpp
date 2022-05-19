#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
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

size_t roundUp(size_t n, int multiple) {
    assert(multiple > 0);
    size_t padded = ((n + multiple - 1) / multiple) * multiple;
    assert(padded >= n);
    assert(padded % multiple == 0);
    return padded;
}

uint32_t* get_n3_counts_v2(size_t N) {
    // Align memory access by padding slightly
    uint32_t* counts = (uint32_t*) aligned_alloc(128, roundUp(sizeof(uint32_t) * (N+1), 128));
    std::fill(counts, counts + N+1, 0);

    // Overflow for uint32 handled by limits on loops + asserts

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
    for (uint32_t i = 1; i <= isqrt(N / 2); i++) {
        uint32_t temp = 2*i*i;
        assert(temp <= N);

        // k = 0
        tuples += 1;
        counts[temp] += 3 * 4;

        // k > 0, k < j
        uint32_t max_k = std::min(i, isqrt(N - temp) + 1);
        for (uint32_t k = 1; k < max_k; k++) {
            uint32_t n = temp + k*k;
            assert(n <= N);
            tuples += 1;
            counts[n] += 24;  // 3 * 8
        }
    }

    // i > j = k
    for (uint32_t i = 1; i <= isqrt(N); i++) {
        uint32_t i_2 = i*i;
        assert(i_2 <= N);

        // j = k = 0
        tuples += 1;
        counts[i_2] += 6;  // 3 * 2

        // j = k, j > 0
        uint32_t max_j = std::min(i, isqrt((N - i_2) / 2) + 1);
        for (uint32_t j = 1; j < max_j; j++) {
            uint32_t n = i_2 + 2*j*j;
            assert(n <= N);
            tuples += 1;
            counts[n] += 24;  // 3 * 8
        }
    }

    for (uint32_t i = 1; i <= isqrt(N); i++) {
        uint32_t i_2 = i*i;
        // i > j, k = 0
        uint32_t max_j = std::min(i, isqrt(N - i_2) + 1);
        for (uint32_t j = 1; j < max_j; j++) {
            uint32_t i_j = i_2 + j*j;
            assert(i_j <= N);
            // k = 0
            tuples += 1;
            counts[i_j] += 24;  // 6 * 4
        }
    }

    // Improve memory access by using smaller counts
    uint16_t* counts_temp = (uint16_t*) aligned_alloc(128, roundUp(sizeof(uint16_t) * (N+1), 128));
    std::fill(counts_temp, counts_temp + N+1, 0);
    uint64_t min_temp = 0;

    uint32_t max_i = isqrt(N);
    // Build sorted list of (j^2 + k^2) to better help with locality
    vector<uint32_t> pairs;
    for (uint32_t i = 3; i <= max_i; i++) {
        uint32_t i_2 = i*i;
        // Largest *VALID* pair (j^2 + k^2)
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
            if ((j_2 + 1) <= max_pair) {
              // Early part of list is constant forever
              {
                  //auto merge_start = std::lower_bound(pairs.begin(), pairs.begin() + start_pairs, j_2 + 1);
                  auto merge_start = std::upper_bound(pairs.begin(), pairs.end(), j_2 + 1);
                  fixed_pos = std::distance(pairs.begin(), merge_start);
                  assert(fixed_pos <= pairs.size());
                  merged = pairs.size() - fixed_pos;

                  /**
                   * For everything up to fixed_pos we can iterate all i's for that item NOW
                   * Do in chunked intervals over both i and pairs so that caching is double good.
                   */
                  /* Have to change max_temp, inplace_merge below when testing this code.
                  size_t pair_end = 0;
                  for (size_t pair_start = 0; pair_start < fixed_pos;) {
                      uint32_t p0 = pairs[pair_start];
                      pair_end = std::distance(
                              pairs.begin(),
                              std::upper_bound(pairs.begin() + pair_start, pairs.begin() + fixed_pos, p0 + 8192));
                      assert(pair_end > pair_start);
                      uint32_t temp_pair_end = pair_end - 1;
                      for (uint32_t later_i = i; later_i <= max_i; later_i++) {
                          uint32_t later_i_2 = later_i * later_i;
                          if (later_i_2 + p0 > N)
                              break;

                          while (later_i_2 + pairs[temp_pair_end] > N) {
                              temp_pair_end -= 1;
                          }

                          for (size_t pi = pair_start; pi <= temp_pair_end; pi++) {
                              uint32_t n = later_i_2 + pairs[pi];
                              assert(n <= N);
                              counts_temp[n] += 1;
                          }
                      }
                      pair_start = pair_end;
                  }
                  pairs.erase(pairs.begin(), pairs.begin() + pair_end);
                  */
              }
              size_t start_pairs = pairs.size();

              uint32_t max_k = std::min(j, isqrt(max_pair - j_2) + 1);
              for (uint32_t k = 1; k < max_k; k++) {
                  uint32_t pair = j_2 + k*k;
                  assert(pair <= max_pair);
                  pairs.push_back(pair);
              }
              assert(pairs.size() > start_pairs);
              pushed = pairs.size() - start_pairs;

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

        assert(pairs.empty() || pairs.back() <= max_pair);

        /**
         * Data being sorted means access to counts is sequential and fast
         *
         * Tried SQL double cached data access (handling batches of pairs & i at same time)
         * Made slow ~10% slower despite several attempts
         *
         * Tried openmp omp parallel with pairs broken into non overlapping intervals, exact same speed.
         * Perf says
         */
        /**
         * Break into wide intervals (guarentee pairs[cpu][-1] < pairs[cpu+1][0]) and use parallel?
         * loop unroll and use sse somehow?
         */
        // if (pairs.size() > 1'000'000) {
        //     /**
        //      * Each interval stop/start should span a multiple of 16.
        //      * Break pairs into groups of 10_000 (plenty of work without to much overhead)
        //      * Threads update [interval[i], interval[i+1])
        //      */
        //     vector<uint32_t> intervals = {0};
        //     size_t start = intervals.back();
        //     while (start < pairs.size()) {
        //         // Start of next interval
        //         start = std::min(start + 10'000, pairs.size());
        //         while (start+1 < pairs.size() && pairs[start]/16 == pairs[start+1]/16) {
        //             start += 1;
        //         }
        //         intervals.push_back(start);
        //     }

        //     #pragma omp parallel for schedule(dynamic)
        //     for (size_t i = 0; i < intervals.size() - 1; i++) {
        //         assert(intervals[i+1] <= pairs.size());
        //         for (uint32_t pi = intervals[i]; pi < intervals[i+1]; pi++) {
        //             uint32_t n = i_2 + pairs[pi];
        //             counts_temp[n] += 1;
        //         }
        //     }
        // } else {
            // Manual loop unrolling of
            // for (auto p : pairs) counts_temp[i_2 + p] += 1;
            uint32_t pi;
            auto counts_start = counts_temp + i_2;
            for (pi = 0; pi < 8 * (pairs.size() / 8); pi += 8) {
                counts_start[pairs[pi]] += 1;
                counts_start[pairs[pi+1]] += 1;
                counts_start[pairs[pi+2]] += 1;
                counts_start[pairs[pi+3]] += 1;
                counts_start[pairs[pi+4]] += 1;
                counts_start[pairs[pi+5]] += 1;
                counts_start[pairs[pi+6]] += 1;
                counts_start[pairs[pi+7]] += 1;
            }
            for (; pi < pairs.size(); pi++) {
                counts_start[pairs[pi]] += 1;
            }
        //}

        /**
         * Numbers may have many representations r3 (https://oeis.org/A025436)
         * Make sure to clear out temp before counts_temp > 65535/4 (or I'll worry about overflow)
         * TODO: could check that sum(counts_temp) = delta tuples
         */
        if (i % 4096 == 0 || i == max_i) {
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

    free(counts_temp);
    printf("\ttuples: %lu\n", tuples);

    return counts;
}

void enumerate_n3(uint64_t N) {
    assert(N <= std::numeric_limits<uint32_t>::max());

    // Would be nice to get incremental results.
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
    vector<int64_t> A000092;
    vector<int64_t> A000223;
    vector<int64_t> A000413;
    for (uint64_t n = 0; n <= N; n++) {
        A_n += counts[n];

        // This can have rounding error at some point.
        double V_n = 4.0/3.0 * M_PI * pow(n, 1.5);
        double P_n = A_n - V_n;
        double record_diff = fabs(P_n) - record;
        if (record_diff > 0) {
            A000092.push_back(n);
            // Can round incorrectly for large n, run validate.py afterwards.
            double P_n_rounded = round(P_n);
            A000223.push_back(P_n_rounded);
            A000413.push_back(A_n);
            record = fabs(P_n);
            uint32_t nth = A000092.size();
            if ((nth < 10) || (nth % 5 == 0) || (nth > 120)) {
                printf("| %3d | %11lu | %14.0f | %15lu | -> %.5f\n", nth, n, P_n_rounded, A_n, record);
            }
        } else if (record_diff > -0.1 || fabs(record_diff) / record < 1e-5) {
            printf("|near | %11lu | miss %9.2f | %15lu | near_record miss by %.5f @ %3lu\n",
                n, P_n, A_n, record - P_n, A000092.size());

        }
    }

    {
        std::ofstream b000092("b000092.txt");
        std::ofstream b000223("b000223.txt");
        std::ofstream b000413("b000413.txt");
        for (size_t i = 0; i < A000092.size(); i++) {
            b000092 << i+1 << " " << A000092[i] << std::endl;
            b000223 << i+1 << " " << A000223[i] << std::endl;
            b000413 << i+1 << " " << A000413[i] << std::endl;
        }
    }
}

int main(void) {
    uint64_t ONE_MILLION = 1'000'000;

    // For 100 terms in 0.12 second
    //enumerate_n3(1560000);

    // For 124 terms in 1.8 seconds
    //enumerate_n3(10 * ONE_MILLION);

    // For 132 terms in 5 seconds
    //enumerate_n3(17 * ONE_MILLION);

    // For 151 terms in 40 seconds
    //enumerate_n3(63 * ONE_MILLION);

    // For 173 terms in 5 minutes
    //enumerate_n3(201 * ONE_MILLION);

    // For 188 terms in 13 minutes
    //enumerate_n3(450 * ONE_MILLION);

    // For 200 terms in 35 minutes
    //enumerate_n3(860 * ONE_MILLION);

    // For 210 terms in 71 minutes
    //enumerate_n3(1400 * ONE_MILLION);

    // For 227 terms in 206 minutes
    //enumerate_n3(2800 * ONE_MILLION);

    // For XXX terms in XXX minutes
    enumerate_n3(4200 * ONE_MILLION);
}
