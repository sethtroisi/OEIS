#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <numeric>
#include <omp.h>
#include <vector>

#define USE_INCR_STRATEGY 1
#define INCREMENT_SIZE (1<<21)

using std::vector;

uint64_t isqrt(uint64_t n) {
    uint64_t t = sqrt(n);
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

/** Optimization breaks on very large function, break edge cases out here */
void handle_small(size_t N, uint32_t* counts, size_t &tuples) {
    // i == j == k
    counts[0] += 1;
    for (uint64_t i = 1; i <= isqrt(N / 3); i++) {
        uint64_t n = 3*i*i;
        assert(n <= N);
        tuples += 1;
        counts[n] += 8;
    }

    // i = j, j > k
    for (uint64_t i = 1; i <= isqrt(N / 2); i++) {
        uint64_t temp = 2*i*i;
        assert(temp <= N);

        // k = 0
        tuples += 1;
        counts[temp] += 3 * 4;

        // k > 0, k < j
        uint64_t max_k = std::min(i-1, isqrt(N - temp));
        for (uint64_t k = 1; k <= max_k; k++) {
            uint64_t n = temp + k*k;
            assert(n <= N);
            tuples += 1;
            counts[n] += 24;  // 3 * 8
        }
    }

    // i > j = k
    for (uint64_t i = 1; i <= isqrt(N); i++) {
        uint64_t i_2 = i*i;
        assert(i_2 <= N);

        // j = k = 0
        tuples += 1;
        counts[i_2] += 6;  // 3 * 2

        // j = k, j > 0
        uint64_t max_j = std::min(i-1, isqrt((N - i_2) / 2));
        for (uint64_t j = 1; j <= max_j; j++) {
            uint64_t n = i_2 + 2*j*j;
            assert(n <= N);
            tuples += 1;
            counts[n] += 24;  // 3 * 8
        }
    }

    for (uint64_t i = 1; i <= isqrt(N); i++) {
        uint64_t i_2 = i*i;
        // i > j, k = 0
        uint64_t max_j = std::min(i-1, isqrt(N - i_2));
        for (uint64_t j = 1; j <= max_j; j++) {
            uint64_t i_j = i_2 + j*j;
            assert(i_j <= N);
            // k = 0
            tuples += 1;
            counts[i_j] += 24;  // 6 * 4
        }
    }
}

void merge_counts(uint32_t* counts, uint16_t* counts_temp, size_t length, bool print_debug) {
    for (uint64_t i = 0; i < length ; i++) {
        counts[i] += (uint32_t) 48 * counts_temp[i];
    }

    uint16_t max_seen = *std::max_element(counts_temp, counts_temp + length);
    std::fill(counts_temp, counts_temp + length, 0);

    if (print_debug || (max_seen > 0x0FFF))
        fprintf(stderr, "\t\tcleared %lu, max_seen=%d\n", length, max_seen);
    assert(max_seen < 0x7FFF); // Make sure we don't overflow
}

uint32_t* get_n3_counts_v2(size_t N) {
    const size_t MAX_PAIR = N * 2 / 3 + 1;

    // Align memory access by padding slightly
    uint32_t* counts = (uint32_t*) aligned_alloc(128, roundUp(sizeof(uint32_t) * (N+1), 128));
    // Improve memory access by using smaller counts
    uint16_t* counts_temp = (uint16_t*) aligned_alloc(128, roundUp(sizeof(uint16_t) * (N+1), 128));
    // Density is ~38%. Slightly better in memory to use pair_count.
    // Readahead would read all of counts anyway and now memory access is perfectly linear.
    uint8_t* pair_count = (uint8_t*) aligned_alloc(128, roundUp(sizeof(uint8_t) * MAX_PAIR, 128));

    std::fill(counts, counts + N+1, 0);
    std::fill(counts_temp, counts_temp + N+1, 0);
    std::fill(pair_count, pair_count + MAX_PAIR, 0);

    /**
     * Pair count eventually overflow uint8, but not within reason. See A025441
     * pair_count[5525] = 6
     * pair_count[2082925] = 18
     * pair_count[243061325] = 48
     * ...
     * pair_count[3,929,086,318,625] = 256
     */

    uint64_t updates_a = 0, updates_b = 0, updates_c = 0;
    uint64_t sum_a, sum_b, sum_c;
    uint64_t min_temp = 0;
    uint64_t num_pairs = 0;

    { // Handle i=j=k, i=j, j=k, k=0, ... cases
        handle_small(N, counts, updates_a);
        sum_a = std::accumulate(counts, counts+N+1, uint64_t{0});
        fprintf(stderr, "\tHandled edge cases, sum: %lu, %lu writes\n\n", sum_a, updates_a);
    }

    { // Build up number of r_2 representations. pair_count[j^2 + k^2] += 1.
        /**
         * Integrate[Min[N - 2*i^2, i^2], {i, 0, Sqrt[N/2]}]
         *     = Integrate[i^2, {i, 0, Sqrt[N/3]}] + Integrate[N - 2*i^2, {i, Sqrt[N/3], Sqrt[N/2]}]
         *     = N^(3/2)/Sqrt[243] + N * (Sqrt[N/2] - Sqrt[N/3]) - 2 * (N^(3/2)/Sqrt[72] - N^(3/2)/Sqrt[243]
         *     = (Sqrt[18] - Sqrt[12]) * N^(3/2) / 9
         */
        float estimated_writes = (sqrt(18) - sqrt(12)) * pow(N, 3.0/2) / 9;

        uint64_t max_i = isqrt(N >> 1);
        assert(2 * (max_i+1)*(max_i+1) > N);
        for (uint64_t i = 3; i <= max_i; i++) {
            uint64_t i_2 = i*i;

            // Largest *VALID* pair (j^2 + k^2)
            // min(N - i_2, (i-1) * (i-1) + (i-2) * (i-2)) = min(N - i*i, 2*i^2 - 6*i + 5)
            const uint64_t max_pair = std::min<uint64_t>(N - i_2, i_2+i_2);
            assert(max_pair < MAX_PAIR);

            uint64_t added = 0;

            // add new pairs for j = i-1
            {
                uint64_t j = i - 1;
                uint64_t j_2 = j*j;
                if ((j_2 + 1) <= max_pair) {

                  uint64_t max_k = std::min(j-1, isqrt(N - i_2 - j_2));
                  assert(j_2 + max_k * max_k <= max_pair);
                  if (max_k < j-1) {
                    assert(j_2 + (max_k+1) * (max_k+1) > max_pair);
                  }
                  auto pair_count_ptr = pair_count + j_2;
                  for (uint64_t k = 1; k <= max_k; k++) {
                      pair_count_ptr[k*k] += 1;
                  }
                  added = max_k;
                  num_pairs += added;
                }
            }

            /**
             * Runs at 1 update / cpu cycle which is probably saturating memory bandwidth
             * One core is doing 4.3e9 updates / sec = 2*2133 MT/s RAM.
             *
             * Weirdly this seems like it might need (read 1B, read 2B, (add), write 2B) / update
             */

            // Only update from the non-constant / changable pairs
            const auto counts_start = counts_temp + i_2;
            for (uint64_t pi = i_2; pi <= max_pair; pi++) {
                counts_start[pi] += pair_count[pi];
            }
            assert(i_2 <= max_pair);
            updates_b += max_pair + 1 - i_2;

            /**
             * Make sure to clear out temp before counts_temp > 65535/4 (or I'll worry about overflow)
             * TODO: could check that sum(counts_temp) = delta tuples
             */
            if (i % 4096 == 0 || i == max_i) {
                assert(min_temp <= i_2);
                uint64_t range = (i_2 - min_temp) + max_pair + 1;
                assert(min_temp + range <= N+1);
                merge_counts(counts + min_temp, counts_temp + min_temp, range, true);
                min_temp = (i+1) * (i+1);

                if (max_pair < i_2) {
                    fprintf(stderr, "\t %5lu/%lu | pairs: %lu (%lu new), %lu writes BREAKING\n",
                            i, max_i, num_pairs, added, updates_b);
                    break;
                }
            }

            if (i && (15 * i) % max_i < 15) {
                uint32_t max_pair_count = *std::max_element(pair_count, pair_count + MAX_PAIR);
                fprintf(stderr, "\t %5lu/%lu | pairs: %lu (max %u, %lu new), %lu writes (%.1f%%)\n",
                        i, max_i, num_pairs, max_pair_count, added,
                        updates_b, 100.0 * updates_b / estimated_writes);
            }
        }

        sum_b = std::accumulate(counts, counts+N+1, uint64_t{0});
        fprintf(stderr, "\tMiddle loop finished, sum: %lu, %lu writes\n\n", sum_b, updates_b);
    }

    assert(std::count(counts_temp, counts_temp + N+1, 0) == (int64_t) (N+1));

    /**
     * pair_count[a] a <= j_2 is constant (this loop changes a >= j_2 + 1)
     *
     * For counts[a] for interval [count_start_point, count_start_point + INCREMENT_SIZE)
     * Handle counts[(i')^2 + pair[...]] += 1 for all i' >= i
     *
     * Because we do multiple i and counts_temp is small should use CPU
     * cache instead of RAM.
     */
    if (!USE_INCR_STRATEGY) {
        uint64_t max_i = isqrt(N);
        for (uint64_t i = 3; i <= max_i; i++) {
            // Add back the constant part of pair_count
            uint64_t i_2 = i*i;
            const auto counts_temp_ptr = counts_temp + i_2;
            uint64_t max_pair = std::min(i_2 - 1, N - i_2);
            //fprintf(stderr, "\t\ti=%lu, [%lu, %lu]\n", i, 0ul, max_pair);
            for (uint64_t pi = 0; pi <= max_pair; pi++) {
                counts_temp_ptr[pi] += pair_count[pi];
            }
            updates_c += max_pair + 1;
            if (i && (15 * i) % max_i < 15) {
                fprintf(stderr, "\t%6lu/%lu pairs: %lu,  %lu writes\n",
                        i, max_i, num_pairs, updates_c);
            }
        }
        merge_counts(counts, counts_temp, N+1, false);
    } else {
        // TODO tune this
        if (omp_get_max_threads() > 4) {
            omp_set_num_threads(4);
        }

        float estimated_writes = (2 - sqrt(2)) * pow(N, 3.0/2) / 3;
        size_t intervals = (N+1-1) / INCREMENT_SIZE + 1;
        #pragma omp parallel for schedule(dynamic, 1)
        for (uint64_t interval = 0; interval < intervals; interval++) {
            // Temp allocate of cache sized interval
            uint16_t thread_counts_temp[INCREMENT_SIZE] = {};

            // Update counts_temp in interval [interval_start, interval_end)
            uint64_t interval_start = interval * INCREMENT_SIZE;
            uint64_t interval_end = std::min(N+1, interval_start + INCREMENT_SIZE);
            uint64_t interval_size = interval_end - interval_start;
            assert(interval_size <= N+1);

            uint64_t min_i = interval_start == 0 ? 3 : isqrt((interval_start >> 1) - 1) + 1;
            uint64_t max_i = isqrt(interval_end-1);

            if (15 * interval % intervals < 15) {
                fprintf(stderr, "\t Interval %lu/%lu [%lu, %lu) with i' [%lu, %lu] %lu writes (%.1f%%)\n",
                        interval, intervals, interval_start, interval_end,
                        min_i, max_i, updates_c, 100.0 * updates_c / estimated_writes);
            }

            if (interval_start > 0) {
                assert(2 * (min_i-1) * (min_i-1) < interval_start);
            }
            assert(2 * min_i * min_i >= interval_start);
            assert(max_i * max_i < interval_end);
            assert((max_i+1) * (max_i+1) >= interval_end);

            // With a large number of i values, chance of overflowing counts_temp.
            const uint64_t I_INCREMENT = 1024;
            for (uint64_t i_range = min_i; i_range <= max_i; i_range += I_INCREMENT) {
                for (uint64_t i = i_range; i < std::min(i_range + I_INCREMENT, max_i + 1); i++) {
                    // Figure out which part of pair_count[pi] is in range
                    uint64_t i_2 = i*i;
                    assert(i_2 < interval_end);

                    uint64_t min_pair = i_2 > interval_start ? 0 : interval_start - i_2;
                    uint64_t max_pair = std::min(i_2, interval_end - i_2);
                    assert(max_pair < MAX_PAIR);
                    assert(min_pair <= max_pair);
                    assert(i_2 + min_pair >= interval_start);
                    assert(i_2 + max_pair-1 < interval_end);
                    // Always use counts_temp [0, interval_size)
                    const auto counts_temp_ptr = thread_counts_temp + i_2 - interval_start;
                    for (uint64_t pi = min_pair; pi < max_pair; pi++) {
                        counts_temp_ptr[pi] += pair_count[pi];
                    }
                    updates_c += max_pair - min_pair;
                }
                merge_counts(counts + interval_start, thread_counts_temp, interval_size, false);
            }
        }
    }

    free(counts_temp);
    free(pair_count);

    sum_c = std::accumulate(counts, counts+N+1, uint64_t{0});
    fprintf(stderr, "\n");
    fprintf(stderr, "\tFinished\n");
    fprintf(stderr, "\twrites: %lu, %lu, %lu\n", updates_a, updates_b, updates_c);
    fprintf(stderr, "\tsums  : %lu, %lu, %lu\n", sum_a, sum_b, sum_c);
    fprintf(stderr, "\n");

    return counts;
}

void enumerate_n3(uint64_t N) {
    /**
     * Memory usage is uint32 + uint16 + uint8 * 2/3 -> 6.66 bytes / N
     * (2^32) * 6.66 -> 27 GB
     */
    assert(N * 20 / 21 <= 48l * (1024*1024*1024));

    /**
     * Everything below i^2 is finalized
     *      Open range is really [i^2, i^2 + (i-1)^2 + (i-2)^2]
     *          Would have to move upper loops into lower loop
     *
     * Pairs has fairly high density
     *      A001481(1000) = 3364 which doesn't account for duplicate representations
     *      Pairs[1000] -> 2664, Pairs[10000] -> 25933
     *      Old appreach stored pairs (4/8 bytes per)
     *      New approach is store 1 byte counter per.
     *
     */

    // Would be nice to get incremental results.
    const auto counts = get_n3_counts_v2(N);

    if (N >= 1000) {
        // Nice verification check from A117609
        uint64_t t = 0;
        for (uint32_t i = 0; i <= 1000; i++) t += counts[i];
        std::cerr << "\tSum of 0..1000: " << t << std::endl;
        assert(t == 132451);

        for (uint64_t i = 1001; i <= N; i++) t += counts[i];
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

        // Double check in case things go wrong
        if ((n % 1'000'000) == 0) {
            fprintf(stderr, "\tA(%lu) = %lu + %u\n", n, A_n, counts[n]);
        }

        // Struggles with rounding error at a(188) and frequently after a(240)
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
            if ((nth < 10) || (nth % 10 == 0) || (nth >= 125)) {
                printf("| %3d | %11lu | %14.0f | %15lu | -> %.5f\n", nth, n, P_n_rounded, A_n, record);
            }
        } else if (record_diff > -0.1 || fabs(record_diff) / record < 1e-5) {
            printf("|near | %11lu | miss %9.2f | %15lu | near_record miss by %.5f @ %3lu\n",
                n, P_n, A_n, record - P_n, A000092.size());

        }
    }
    free(counts);

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

    // For 100 terms in 0.07 second
    //enumerate_n3(1560000);

    // For 124 terms in 0.8 seconds
    //enumerate_n3(10 * ONE_MILLION);

    // For 129 terms in 1.7 seconds
    //enumerate_n3(17 * ONE_MILLION);

    // For 138 terms in 8 seconds
    //enumerate_n3(40 * ONE_MILLION);

    // For 151 terms in 16 seconds
    //enumerate_n3(63 * ONE_MILLION);

    // For 158 terms in 32 seconds
    //enumerate_n3(100 * ONE_MILLION);

    // For 170 terms in 90 seconds
    //enumerate_n3(201 * ONE_MILLION);

    // For 188 terms in 5 minutes
    //enumerate_n3(450 * ONE_MILLION);

    // For 200 terms in 13 minutes
    //enumerate_n3(860 * ONE_MILLION);

    // For 210 terms in 26 minutes
    //enumerate_n3(1400 * ONE_MILLION);

    // For 227 terms in 78 minutes
    //enumerate_n3(2800 * ONE_MILLION);

    // For 235 terms in 133 minutes
    //enumerate_n3(4200 * ONE_MILLION);

    // For 240 terms in 167 minutes
    //enumerate_n3(4800 * ONE_MILLION);

    // For 253 terms in ~230 minutes
    //enumerate_n3(6000 * ONE_MILLION);

    // No terms from 5967m to 6500m
    // Result was off by 3 * 256 so possible an overflow in pair_count?
    enumerate_n3(6000 * ONE_MILLION);

    return ONE_MILLION - ONE_MILLION;
}
