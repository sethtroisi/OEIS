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

void merge_counts(uint32_t* counts, uint16_t* counts_temp, size_t length) {
    for (uint64_t i = 0; i < length ; i++) {
        counts[i] += (uint32_t) 48 * counts_temp[i];
    }

    uint16_t max_seen = *std::max_element(counts_temp, counts_temp + length);
    std::fill(counts_temp, counts_temp + length, 0);

    fprintf(stderr, "\t\tcleared %lu, max_seen=%d\n", length, max_seen);
    assert(max_seen < 0x3FFF); // Make sure we don't overflow
}

uint32_t* get_n3_counts_v2(size_t N) {
    // Align memory access by padding slightly

    uint32_t* counts = (uint32_t*) aligned_alloc(128, roundUp(sizeof(uint32_t) * (N+1), 128));
    // Improve memory access by using smaller counts
    uint16_t* counts_temp = (uint16_t*) aligned_alloc(128, roundUp(sizeof(uint16_t) * (N+1), 128));
    // Density is ~38%. Slightly better in memory to use pair_count.
    // Readahead would read all of counts anyway and now memory access is perfectly linear.
    uint8_t* pair_count = (uint8_t*) aligned_alloc(128, roundUp(sizeof(uint8_t) * (N+1), 128));

    std::fill(counts, counts + N+1, 0);
    std::fill(counts_temp, counts_temp + N+1, 0);
    std::fill(pair_count, pair_count + N+1, 0);

    uint64_t updates = 0;
    uint64_t min_temp = 0;
    uint64_t num_pairs = 0;

    handle_small(N, counts, updates);

    uint64_t max_i = isqrt(N);
    for (uint64_t i = 3; i <= max_i; i++) {
        uint64_t i_2 = i*i;

        // Largest *VALID* pair (j^2 + k^2)
        // min(N - i_2, (i-1) * (i-1) + (i-2) * (i-2)) = min(N - i*i, 2*i^2 - 6*i + 5
        const uint64_t max_pair = std::min<uint64_t>(N - i_2, i_2+i_2); //(i-1)*(i-1) + (i-2)*(i-2));

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
        updates += i_2 > max_pair ? 0 : max_pair + 1 - i_2;

        /**
         * Make sure to clear out temp before counts_temp > 65535/4 (or I'll worry about overflow)
         * TODO: could check that sum(counts_temp) = delta tuples
         */
        if (i % 4096 == 0 || i == max_i) {
            assert(min_temp <= i_2);
            uint64_t range = (i_2 - min_temp) + max_pair + 1;
            assert(min_temp + range <= N+1);
            merge_counts(counts + min_temp, counts_temp + min_temp, range);
            min_temp = (i+1) * (i+1);
        }

        if (i && (20 * i) % max_i < 20) {
            fprintf(stderr, "\t%6lu/%lu pairs: %lu (%lu new), %lu writes\n",
                    i, max_i, num_pairs, added, updates);
        }
    }

    fprintf(stderr, "\n");

    /**
     * pair_count[a] a <= j_2 is constant (this loop changes a >= j_2 + 1)
     *
     * For counts[a] for interval [count_start_point, count_start_point + INCREMENT_SIZE)
     * Handle counts[(i')^2 + pair[...]] += 1 for all i' >= i
     *
     * Because we do multiple i and counts_temp is small should use CPU
     * cache instead of RAM.
     */
    for (uint64_t i = 3; i <= max_i; i++) {
        // Add back the constant part of pair_count
        uint64_t i_2 = i*i;
        const auto counts_temp_ptr = counts_temp + i_2;
        uint64_t max_pair = std::min(i_2-1, N - i_2);
        for (uint64_t pi = 0; pi <= max_pair; pi++) {
            counts_temp_ptr[pi] += pair_count[pi];
        }
        updates += max_pair + 1;

        if (i && (20 * i) % max_i < 20) {
            fprintf(stderr, "\t%6lu/%lu pairs: %lu,  %lu writes\n",
                    i, max_i, num_pairs, updates);
         }
    }
    merge_counts(counts, counts_temp, N);


    free(counts_temp);

    return counts;
}

void enumerate_n3(uint64_t N) {
    /**
     * Memory usage is uint32 + uint16 + uint8 -> 7 bytes / N
     * (2^32) * 7 -> 28 GB
     */
    assert(N * 7 <= 32l * (1024*1024*1024));

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

    if (N > 1000) {
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
            if ((nth < 10) || (nth % 10 == 0) || (nth > 150)) {
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

    // For 129 terms in 2.4 seconds
    //enumerate_n3(17 * ONE_MILLION);

    // For 138 terms in 15 seconds
    enumerate_n3(40 * ONE_MILLION);

    // For 151 terms in 32 seconds
    //enumerate_n3(63 * ONE_MILLION);

    // For 158 terms in 62 seconds
    //enumerate_n3(100 * ONE_MILLION);

    // For 170 terms in 3 minutes
    //enumerate_n3(201 * ONE_MILLION);

    // For 188 terms in 11 minutes
    //enumerate_n3(450 * ONE_MILLION);

    // For 200 terms in 27 minutes
    //enumerate_n3(860 * ONE_MILLION);

    // For 210 terms in 54 minutes
    //enumerate_n3(1400 * ONE_MILLION);

    // For 227 terms in <206 minutes
    //enumerate_n3(2800 * ONE_MILLION);

    // For 235 terms in 290 minutes
    //enumerate_n3(4200 * ONE_MILLION);

    return ONE_MILLION - ONE_MILLION;
}
