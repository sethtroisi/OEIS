// g++ -g -O3 --std=c++17 -Werror -Wall A217259.cpp -lgmpxx -lgmp -pthread -fopenmp && time ./a.out

/*
11 threads

225569500 	100,580,280,191 		150.0 seconds elapsed 670.5M/s
434699100 	205,449,033,509 		340.0 seconds elapsed 604.3M/s
625725600 	305,209,486,511 		540.0 seconds elapsed 565.2M/s
994193100 	504,381,603,587 		980.0 seconds elapsed 514.7M/s
1875582200	1,002,891,288,197		2260.0 seconds elapsed 443.8M/s
3556689700	2,002,379,757,239		5400.0 seconds elapsed 370.8M/s
8316959500	5,002,892,027,909		18100.0 seconds elapsed 276.4M/s
15838597700	10,002,668,115,959		46320.0 seconds elapsed 215.9M/s
30198403100	19,999,673,308,007		122420.0 seconds elapsed 163.4M/s

twin prime count (up to 2e13) confirmed by primesieve in 432 seconds!

TODO
Add 2541865828322 to A063680 and update limit when done
*/

#include <algorithm>
#include <cassert>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <memory>
#include <mutex>
#include <thread>
#include <utility>
#include <vector>

#include <gmp.h>
#include <gmpxx.h>

using std::vector;

uint64_t calc_isqrt(uint64_t n) {
    // isqrt with gmp, could use faster uint64_t math
    mpz_class sqrt = n;
    mpz_sqrt(sqrt.get_mpz_t(), sqrt.get_mpz_t());
    return mpz_get_ui(sqrt.get_mpz_t());
}

/**
 * Calculate sum of divisors (excluding 1 and self) for numbers in [start, start+n)
 *
 * Ideas on how to improve
 *
 * Return  low bits (uint16_t/uint32_t) of sigma
 *      Takes less than O(1M) to test false positives
 *
 * Only look at primes (10x fewer "factors")
 *      multiply (1 + p1 + p1^2 + p2^3) * (1 + p2 + p2^3) * (1 + p3)
 *      In theory a wheel like (2^5 3^5 5^2) makes this easier
 *
 *      Notes from trying this:
 *          Everything was hard and involved division.
 *          Tried `count_factor(index, prime)` and also keeping an array of `[n % mod_p^i]`
 *              But both are slow (in different ways)
 *          Tried doing multiple passes
 *              Never finished this code but looked like it was going to involve the same amount of division
 *                  Would have to divide off previous mult for an index (back to needing count_factor[index, prime])
 *          Thought about in single-threaded variant that tracks [start % mod_p^i] but [p * powers]
 *              seems like I would be keeping track of the same number of elements
 *                  O(n/log(n) * log(n)) vs than just the summation strategy of O(7 * n) (7 = sum(1/i, i=15...sqrt(n))) .
 *                      seems like eventually would be fewer updates?
 *          Worried that I'm not covering large primes
 *              decided I would either need to track ALL primes up to `start` e.g 10^9+
 *              OR would need to keep division of `remainder[i] = (start + i) / prime` for primes up to sqrt
 *                  then would have to multiply counts[i] * (1 + remainder[i])
 *
 *          Everything seems division heavy and current approach is division light and very fast.
 *              I think I'm producing a sigma every 20-80 cycles based on 4Ghz / 130M/s = 31cycles / i
 */
vector<uint64_t> SegmentedSieveOfSigma(uint64_t start, uint64_t N) {
    auto past = start + N;
    auto sums = vector<uint64_t>(N, 0);

    /*
    // Adjust to include n & 1 as a divisor of n
    for (uint64_t i = 0; i < N; i++) {
        sums[i] += start + i + 1;
    }
    */

    if (start == 0) {
        sums[0] = -3;  // undefined
        sums[1] = -5;  // different undefined
    }

    uint64_t isqrt = calc_isqrt(past - 1);
    assert( isqrt * isqrt < past );
    assert( (isqrt+1) * (isqrt+1) >= past );

    /**
     * Handle the few factors with start <= f^2 < start+N seperatly
     */
    for (; isqrt >= 2; isqrt--) {
        uint32_t factor = isqrt;
        uint64_t f2 = (uint64_t) factor * factor;
        if (f2 < start)
          break;

        assert(f2 < start + N);

        uint32_t index = f2 - start;
        // n = factor^2 only adds factor because n / factor = factor.
        sums[index] += factor;

        // n = factor * (factor+1)
        uint32_t aliquot_add = (factor) + (factor + 1);
        for (index += factor; index < N; index += factor) {
            sums[index] += aliquot_add++;
        }
    }

    uint64_t start_m_1 = start - 1;

    // Hand unrolled loops for very small factor
    uint64_t factor = 2;
    for (; factor <= std::min(isqrt, N/5); factor++) {
        uint64_t count = start_m_1 / factor + 1;
        uint32_t index = count * factor - start;
        uint64_t add = factor + count;

        // Loop unrolled 4x for small factors
        for (; index + (factor<<2) < N; ) {
            sums[index           ]  += add++;
            sums[index +   factor]  += add++;
            sums[index + 2*factor]  += add++;
            sums[index + 3*factor]  += add++;
            index += factor<<2;
        }

        for (; index < N; index += factor)
            sums[index] += add++;
    }

    // Handles factors that can appear more than once
    for (; factor <= std::min(isqrt, N); factor++) {
        uint64_t count = start_m_1 / factor + 1;
        uint32_t index = count * factor - start;
        uint64_t add = factor + count;

        for (; index < N; index += factor)
            sums[index] += add++;
    }

    // Handles larger factors which can only appear once
    for (; factor <= isqrt; factor++) {
        uint64_t count = start_m_1 / factor + 1;
        uint32_t index = count * factor - start;
        uint64_t add = factor + count;

        if (index < N)
          sums[index] += add;
    }

    return sums;
}

class SegmentedSieveSigma {
    public:
        SegmentedSieveSigma(uint64_t start, uint64_t length) : sieve_length(length) {
            assert(sieve_length == 2*3*2*5*1*7*2*3*1*11*1*13*1*1);
            assert(sieve_length == 360360);

            base_counts.resize(sieve_length);
            base_delta.resize(sieve_length);

            // So that we don't have to deal with f^2 entries after 0th interval.
            assert(MAX_BASE_F*MAX_BASE_F < sieve_length);

            for (uint32_t f = 2; f <= MAX_BASE_F; f++) {
                if (sieve_length % f == 0) {
                    // Account for base_delta being added once
                    for (uint32_t i = 0; i < sieve_length; i += f)
                        base_counts[i] += f + (i / f);

                    for (uint32_t i = 0; i < sieve_length; i += f)
                        base_delta[i] += sieve_length / f;
                }
            }

            sums.resize(sieve_length);
            if (start > 0) {
                jump_to(start);
            }
        }

        void jump_to(uint64_t new_start);
        const vector<uint64_t>& next(uint64_t verify_start);

    private:
        const uint32_t sieve_length;

        // Results are stored here!
        vector<uint64_t> sums;

        // Start of the next interval, must be a multiple of sieve_length
        uint64_t start = 0;

        // counts[i] = ceil(start / i)
        vector<uint64_t> counts = {0,0};

        /**
         * Handle very small factors with a wheel approach
         * wheel size = lcm(2 ... 15) = 360360
         *
         * The increment (base_delta) in each position is sum(d / wheel_size, d | 360360)
         * Not that we limit d to 600, to avoid complex handling of d^2
         *
         * sum( 1/i, i >= 2 && i | 360360 ) = 3.1
         *
         *      Means reduction of 3.1 * sieve_length random writes to 2*N linear writes
         *          base_counts[i] += base_delta[i]
         *          sum[i] = base_counts[i]
         *
         * sum( 1/i, 2 < i < sieve_length ) = 13
         *      Reducing 1/4 of that is nice!
         */
        vector<uint64_t> base_counts;
        vector<uint64_t> base_delta;
        const uint32_t MAX_BASE_F = 600; // isqrt(360360-1)
};


void SegmentedSieveSigma::jump_to(uint64_t new_start) {
    assert(new_start > start);
    assert((new_start - start) % sieve_length == 0);
    assert(new_start % sieve_length == 0);
    uint64_t segment_jump = (new_start - start) / sieve_length;
    start = new_start;

    // Add factors, f, while f*f < start
    auto isqrt = calc_isqrt(start - 1);
    assert((isqrt+1) >= counts.size());

    // Possible add new empty counts
    counts.resize(isqrt + 1);

    for (uint64_t f = 2; f <= isqrt; f++) {
        // -start % f
        uint64_t count = (start-1) / f + 1;

        if (f <= MAX_BASE_F && sieve_length % f == 0)
            counts[f] = -1ul; // should break code if encountered
        else
            counts[f] = count;
    }

    for (uint64_t i = 0; i < sieve_length; i++) {
        base_counts[i] += segment_jump * base_delta[i];
    }

    // Spot check
    assert(new_start > 0);
    assert(base_counts[3] == 3 + (new_start + 3) / 3);
    assert(base_counts[7] == 7 + (new_start + 7) / 7);
    assert(base_counts[13] == 13 + (new_start + 13) / 13);
 }

const vector<uint64_t>& SegmentedSieveSigma::next(uint64_t verify_start) {
    assert(start == verify_start);

    // First interval is WEIRD with many f^2 < sieve_length
    if (start == 0) {
        auto t = SegmentedSieveOfSigma(0, sieve_length);
        for (uint64_t i = 0; i < sieve_length; i++) sums[i] = t[i];

        // let jump_to() correctly update state
        jump_to(sieve_length);
        return sums;
    }

    //auto sums = vector<uint64_t>(sieve_length, 0);
    //std::memset(&sums[0], 0, sieve_length * sizeof(sums[0]));

    for (uint64_t i = 0; i < sieve_length; i++) {
        sums[i] = base_counts[i];
        base_counts[i] += base_delta[i];
    }

    auto past = start + sieve_length;

    // For all new factors
    for(uint64_t f = counts.size(), f2 = f * f; f2 < past; f++, f2 = f*f) {
        assert(f2 >= start);

        uint32_t index = f2 - start;
        assert(index < sieve_length);

        // n = factor^2 only adds factor because n / factor = factor.
        sums[index] += f;

        // Add offset for lower loop to handle
        assert(counts.size() == f);

        if (f <= MAX_BASE_F && sieve_length % f == 0)
            counts.push_back(-1ul); // should break code if encountered
        else
            counts.push_back(f + 1);
    }
    assert((uint64_t) (counts.size()+1)*(counts.size()+1) >= past);

    // Hand unrolled loops for very small factor
    uint64_t factor = 2;
    for (; factor < std::min<uint32_t>(sieve_length / 9, counts.size()); factor++) {
        if (factor <= MAX_BASE_F && sieve_length % factor == 0)
            continue;

        auto count = counts[factor];
        assert(count != -1ul);
        //assert(count  == std::max(factor + 1, ((start + factor - 1) / factor)));
        uint32_t index = count * factor - start;
        auto add = count + factor;

        if (index >= sieve_length) {
            // can happen when f is first added to the array (mostly in first interval)
            continue;
        }

        // ceil((sieve_length - index) / factor)
        int32_t updates = (sieve_length - index - 1) / factor + 1;
        assert(index + (updates-1) * factor < sieve_length);
        assert(index + updates * factor >= sieve_length);

        // Loop unrolled 8x for small factors
        for (; updates >= 8; updates -= 8) {
            sums[index           ]  += add++;
            sums[index +   factor]  += add++;
            sums[index + 2*factor]  += add++;
            sums[index + 3*factor]  += add++;
            sums[index + 4*factor]  += add++;
            sums[index + 5*factor]  += add++;
            sums[index + 6*factor]  += add++;
            sums[index + 7*factor]  += add++;
            index += factor<<3;
        }

        for (; updates > 0; updates--) {
            sums[index] += add++;
            index += factor;
        }

        assert((uint32_t) index >= sieve_length);
        counts[factor] = add - factor;
    }

    for (; factor < std::min<uint32_t>(sieve_length / 2, counts.size()); factor++) {
        auto count = counts[factor];
        assert(count != -1ul);
        //assert(count  == std::max(factor + 1, ((start + factor - 1) / factor)));
        uint32_t index = count * factor - start;
        auto add = count + factor;

        for (; index < sieve_length; index += factor)
            sums[index] += add++;

        counts[factor] = add - factor;
    }

    // Handles factors that appear at least once (but possible more times)
    for (; factor < std::min<uint32_t>(counts.size(), sieve_length); factor++) {
        auto count = counts[factor];
        assert(count != -1ul);
        uint32_t index = count * factor - start;
        auto add = count + factor;

        for (; index < sieve_length; index += factor)
            sums[index] += add++;

        counts[factor] = add - factor;
    }

    // Handles larger factors that can only appear once
    for (; factor < counts.size(); factor++) {
        auto count = counts[factor];
        assert(count != -1ul);
        uint32_t index = count * factor - start;

        if (index < sieve_length) {
            sums[index] += count + factor;
            counts[factor]++;
        }
    }

    start += sieve_length;
    return sums;
}

bool sigmaSelfCheck() {
    uint32_t SEGMENT = 360360;

    auto errors = 0;
    uint64_t summation = 0;

    auto sieve = SegmentedSieveSigma(0, SEGMENT);

    for (uint32_t start = 0; start < 3*SEGMENT; start += SEGMENT) {
        auto sigmas_a = sieve.next(start);
        auto sigmas_b = SegmentedSieveOfSigma(start, SEGMENT);

        assert(sigmas_a.size() == SEGMENT);
        assert(sigmas_b.size() == SEGMENT);

        uint32_t first_i = (start > 0) ? 0 : 2;

        for (uint32_t i = first_i; i < SEGMENT; i++) {
            summation += sigmas_a[i] + (start + i + 1);

            if (sigmas_a[i] != sigmas_b[i]) {
                printf("Self Check mismatch at %u (%u) | %lu vs %lu\n",
                        start + i, i, sigmas_a[i], sigmas_b[i]);
                errors += 1;
                if (errors > 5)
                    return false;
            }
        }
        //printf("\tSelf check summation [0, %u) = %lu\n", start + SEGMENT, summation);
    }

    // 106804038823, 427218070633, 961242021633
    uint64_t expected = 961242021633;
    if (summation != expected) {
        printf("\tSummation %lu != %lu !!!\n", summation, expected);
        assert(summation == expected);
    }
    //printf("\n");

    return errors == 0;
}

class A217259 {
    public:
        A217259(uint64_t start, uint64_t stop, uint64_t n)
            : START(start), STOP(stop), SEGMENT(n) {}

        bool test_composite_match(int16_t dist, uint64_t i);
        void print_match(int16_t dist, uint64_t i);
        void print_match_and_test(int16_t dist, uint64_t i) {
            print_match(dist, i);
            if (dist > 0)
                assert(test_composite_match(dist, i));
        }
        void print_results(uint64_t last_included);

        void iterate();
        void multithreaded_iterate() {
            std::thread t1(&A217259::worker_coordinator, this);
            std::thread t2(&A217259::worker_thread, this);

            t1.join();
            t2.join();
        }

    private:
        void worker_thread();
        void worker_coordinator();

        inline void process_pair(uint32_t d, uint64_t a, uint64_t sigma_a, uint64_t sigma_b);

        const uint64_t START;
        const uint64_t STOP;
        const uint64_t SEGMENT;
        // DIST = 2 (A217259), 6 (A054903), 7 (A063680), 8 (A059118), 12 (A054902)
        static const uint32_t MAX_DIST = 11;

        uint64_t last_match = 0;
        int64_t found_prime[MAX_DIST+1] = {};
        int64_t total_composites = 0;
        int64_t found_composite[MAX_DIST+1] = {};
        int64_t print_mult = 1;
        uint64_t next_time = 5;

        std::chrono::time_point<std::chrono::system_clock> S =
            std::chrono::system_clock::now();

        // Guard for results
        std::mutex g_control;
        std::condition_variable work_ready;

        // Multithreaded results output
        typedef vector<std::pair<int8_t,uint64_t>> terms_t;
        typedef std::unique_ptr<terms_t> terms_ptr_t;
        vector<std::pair<uint64_t,terms_ptr_t>> results;
};


void A217259::print_results(uint64_t last_stop) {
    std::chrono::duration<double> elapsed = std::chrono::system_clock::now() - S;

    printf("\n\n\n\n");
    printf("Results From [%lu, %lu] in groups of %lu\ttook %.0f seconds",
            START, last_stop, SEGMENT, elapsed.count());

    // TODO track if we ran singlethreaded or multithreaded
    printf("\n");
    printf("Found %'lu primes.\n", found_prime[0]);
    printf("Found prime pairs at these distances:\n");
    for (uint32_t d = 1; d <= MAX_DIST; d++) {
        if (found_prime[d] > 1)
            printf("\t%d -> %'lu pairs\n", d, found_prime[d]);
    }
    printf("Found %lu total composites, from these distances:\n", total_composites);
    for (uint32_t d = 1; d <= MAX_DIST; d++) {
        if (found_composite[d])
            printf("\t%d -> %'lu pairs\n", d, found_composite[d]);
    }
}

bool A217259::test_composite_match(int16_t dist, uint64_t i) {
    // Already known terms
    // DIST = 2 (A050507), 6 (A054903), 7 (A063680), 8 (A059118), 12 (A054902)
    if (dist == 2)
        for (uint64_t known : {434, 8575, 8825})
            if (i == known)
                return true;

    if (dist == 4)
        for (uint64_t known : {
                305635357,})
            if (i == known)
                return true;

    if (dist == 6)
        for (uint64_t known : std::initializer_list<uint64_t>{
                104, 147, 596, 1415, 4850, 5337, 370047, // 1630622,
                35020303, 120221396, 3954451796, 742514284703})
            if (i == known)
                return true;

    if (dist == 7)
        for (uint64_t known : std::initializer_list<uint64_t>{
                74, 531434, 387420482
                })
            if (i == known)
                return true;

    if (dist == 8)
        for (uint64_t known : std::initializer_list<uint64_t>{
                27,1615,1885,218984,4218475,312016315,746314601,
                1125845307,1132343549,1296114929,9016730984,
                303419868239,1197056419121,2065971192041,
                2948269852109,4562970154601
                })
            if (i == known)
                return true;


    if (dist == 10)
        for (uint64_t known : std::initializer_list<uint64_t>{
                195556, 1152136225
                })
            if (i == known)
                return true;


    if (dist == 12)
        for (uint64_t known : std::initializer_list<uint64_t>{
                65,170,209,1394,3393,4407,4556,11009,13736,27674,
                38009,38845,47402,76994,157994,162393,184740,
                186686,209294,680609,825359,954521,1243574,
                2205209,3515609,4347209,5968502,6539102,6916241,
                8165294,10352294,10595009,10786814
                })
            if (i == known)
                return true;

    for (uint64_t test : {i, i + dist}) {
        mpz_class t = test;
        int primality = mpz_probab_prime_p(t.get_mpz_t(), 20);
        if (primality != 0) {
            printf("\n\nPRIMALITY MISMATCH FOR %lu | %d vs %d\n\n", test, primality, dist);
            exit(1);
        }

        if (primality == 0 && test == i) {
            printf("\n\n\n\nNEW! %u, %lu | %lu is composite(%d)!\n\n\n\n\n",
                    abs(dist), i, mpz_get_ui(t.get_mpz_t()), primality);
        }
    }
    return true;
}

void A217259::print_match(int16_t dist, uint64_t i) {
    // Don't duplicate matches like [3,5],  [3,7], [3,11]

    // A composite!
    if (dist > 0) {
        total_composites += 1;
        found_composite[dist] += 1;
    } else {
        found_prime[-dist] += 1;

        if (dist == -1)
            last_match = i;
    }

    if (dist > 0) {
        printf("%2ld composite\t%'-16lu \t\t%ld'th composite with d=%d\n",
                total_composites, i, found_composite[dist], dist);
    } else if (dist == -2) {
        auto twin_primes = found_prime[2];

        // So we can verify against A146214
        if (twin_primes == 6*print_mult)
            print_mult *= 10;

        if (twin_primes % print_mult == 0) {
            printf("%-10ld\t%'-16lu \t\ttwin prime\n", twin_primes, i);
        } else if (twin_primes % 10000 == 0) {
            // Avoid calling sys_clock on every find.
            std::chrono::duration<double> elapsed = std::chrono::system_clock::now() - S;
            if (elapsed.count() > next_time) {
                // TODO add last interval rate
                float rate = (i - START) / elapsed.count() / 1e6;
                printf("%-10ld\t%'-16lu\t\t%.1f seconds elapsed %.1fM/s\n",
                        twin_primes, i, elapsed.count(), rate);

                // 5,10,15 ... 95,100,110,120,130 ... 190,200,220,240 ... 960,1000,1080,1160
                next_time += 5 * (1 + (next_time >= 30)) * (1 + (next_time >= 200)) * (1 + (next_time >= 1000));
            }
        }
    }
}

inline void A217259::process_pair(
        uint32_t d, uint64_t a, uint64_t sigma_a, uint64_t sigma_b) {

    // When sigma is sigma_1 (includes n and 1)
    //if (sigma_a + d == sigma_b) {
    //    if (sigma_a == 0 && sigma_b == 0) {
    //        print_match(-d, a);

    // When sigma is (sigma_1 - 1 - n) == 0 for primes
    if (sigma_a == sigma_b) {
        if (a > STOP) {
            printf("Found match(%lu with %lu) just past STOP\n", sigma_a, sigma_b);
            return;
        }

        if (sigma_a == 0 && sigma_b == 0) {
            print_match(-d, a);
        } else {
            printf("\n\tsigmas[%lu] = %lu, sigmas[%lu] = %lu\n",
                a, sigma_a, a + d, sigma_b);
            print_match_and_test(d, a);
        }
    }
}

void A217259::iterate() {
    // EXCLUDES STOP

    uint64_t last_sigmas[MAX_DIST];
    for (size_t i = 0; i < MAX_DIST; i++) {
        last_sigmas[i] = -200ul;
    }

    auto sieve = SegmentedSieveSigma(START, SEGMENT);

    uint64_t start;
    for (start = START; start < (STOP + MAX_DIST); start += SEGMENT) {
        auto sigmas2 = SegmentedSieveOfSigma(start, SEGMENT);
        auto sigmas = sieve.next(start);

        for (uint32_t i = 0; i < SEGMENT; i++)
            if (sigmas[i] != sigmas2[i])
                printf("\n\nSIGMA MISMATCH @ %lu (%u) | %lu vs %lu\n\n\n",
                        start + i, i, sigmas[i], sigmas2[i]);

        // Considering keep sum of sigma as a spot check

        // test i less than or equal to these values in the three loops
        uint64_t last_i_head = MAX_DIST-1;
        uint64_t last_i_main = SEGMENT-MAX_DIST-1;
        uint64_t last_i_tail = SEGMENT-1;

        if (start + SEGMENT <= STOP) {
            // Run the full interval!
        } else if (STOP <= start) {
            // Do some portion of head search, none of main body
            last_i_head = (MAX_DIST-1) - (start - STOP);
            last_i_main = 0;
            last_i_tail = 0;
            assert(start - MAX_DIST + last_i_head == STOP - 1);
            assert(last_i_head <= MAX_DIST-1);

        } else if (STOP < (start + SEGMENT)) {
            // Do some portion of the main search

            auto last_i = STOP - start - 1;
            assert(last_i >= 0);
            assert(last_i < SEGMENT-1);

            // Can do some portion of main search
            if (last_i < last_i_main) {
                last_i_main = std::min(last_i, last_i_main);
                assert(start + last_i_main == (STOP - 1));
            }
            last_i_tail = last_i;

            assert(last_i_main + MAX_DIST < SEGMENT);
            assert(start + last_i_tail == (STOP-1));

        } else {
            printf("??? ERROR: start: %lu, segment: %lu, stop: %lu, max_dist: %u\n",
                    start, SEGMENT, STOP, MAX_DIST);
            exit(2);
        }

        assert(start == 0 || (start - MAX_DIST + last_i_head) < STOP);
        assert(last_i_main == 0 || last_i_main + MAX_DIST < SEGMENT);
        assert(last_i_main == 0 || last_i_main + start < STOP);
        assert(last_i_tail == 0 || last_i_tail + start < STOP);

        if (start > 0) {
            for (uint32_t i = 0; i <= last_i_head; i++) {
                for (uint32_t d = MAX_DIST - i; d <= MAX_DIST; d++) {
                    // 0 is start-MAX_DIST   | which is MAX_DIST distance from 0
                    // 1 is start-MAX_DIST+1 | which is MAX_DIST-1 distance from 0
                    // 1 is start-MAX_DIST+1 | which is MAX_DIST   distance from 1
                    // 2 is start-MAX_DIST+2 | which is MAX_DIST-2 distance from 0
                    // 2 is start-MAX_DIST+2 | which is MAX_DIST-1 distance from 1
                    // 2 is start-MAX_DIST+2 | which is MAX_DIST   distance from 2
                    //assert(i + d >= MAX_DIST);

                    process_pair(d, start - MAX_DIST + i, last_sigmas[i], sigmas[i + d - MAX_DIST]);
                }
            }
        }

        // MAIN checking
        for (uint32_t i = 0; i <= last_i_main; i++) {
            if(sigmas[i] == 0)
                found_prime[0] += 1;

            for (uint32_t d = 2; d <= MAX_DIST; d++) {
                process_pair(d, start + i, sigmas[i], sigmas[i + d]);
            }
        }

        // Checking of final items which may roll over to next window.
        for (uint32_t i = SEGMENT - MAX_DIST; i <= last_i_tail; i++) {
            last_sigmas[i - (SEGMENT-MAX_DIST)] = sigmas[i];

            if(sigmas[i] == 0)
                found_prime[0] += 1;

            // All the pairs which stay inside [start, start + SEGMENT)
            for (uint32_t d = 2; d < (SEGMENT - i); d++) {
                process_pair(d, start + i, sigmas[i], sigmas[i + d]);
            }
        }
    }

    print_results(STOP-1);
}


void A217259::worker_thread() {
    #pragma omp parallel for schedule(dynamic, 4)
    for (uint64_t start = START; start < STOP; start += (SEGMENT - MAX_DIST)) {
        // Calculate results
        vector<uint64_t> sigmas = SegmentedSieveOfSigma(start, SEGMENT);

        terms_ptr_t terms(new terms_t());
        uint32_t found_primes = 0;

        auto stop_i = SEGMENT-MAX_DIST;
        if (start + stop_i > STOP) {
            assert(start < STOP);
            stop_i = (STOP - start);
            assert(stop_i + MAX_DIST <= SEGMENT);
        }
        //printf("\t[%lu, %lu)\n", start, start + stop_i);

        // Only look at i in the non-overlapping part
        for (uint32_t i = 0; i < stop_i; i++) {
            if (sigmas[i] == 0)
                found_primes += 1;

            for (uint16_t dist = 2; dist <= MAX_DIST; dist++) {
                if (sigmas[i] == sigmas[i + dist]) {
                    /**
                     * Old code verified i+1 and i-1 are prime explicitly.
                     * It's instant to instead check sigmas[i] == sigmas[j] == 0
                     * Can spot check with twin prime count (A146214)
                     */
                    if (sigmas[i] == 0 && sigmas[i+dist] == 0) {
                        // Prime pair (e.g. twin, cousin, ...)
                        terms->push_back({-dist, start + i});
                    } else {
                        printf("\n\tsigmas[%lu] = %lu, sigmas[%lu] = %lu\n",
                            start+i, sigmas[i], start+i+dist, sigmas[i+dist]);

                        assert(test_composite_match(dist, start + i));
                        terms->push_back({dist, start + i});
                    }
                }
            }
        }

        // Wait till I can safely queue results.
        std::unique_lock<std::mutex> guard(g_control);

        //printf("\t\tWork %lu -> %lu ready\n", start, terms.size());
        results.emplace_back(start, std::move(terms));
        found_prime[0] += found_primes;

        // Let iteratator advance
        guard.unlock();
        work_ready.notify_one();
    }
    printf("\t\tAll work finished!\n");
}

void A217259::worker_coordinator() {
    // Wait for added item
    std::unique_lock<std::mutex> guard(g_control);
    // Wait for some work ready
    work_ready.wait(guard);

    uint64_t start = START;
    while (start < STOP) {
        if (!guard) {
            guard.lock();
        }

        if (results.empty()) {
            // Release lock and wait for worker_thread
            work_ready.wait(guard);
            continue;
        }

        terms_ptr_t term_ptr;

        // Check if start in results
        for (size_t i = 0; i < results.size(); i++) {
            if (results[i].first == start) {
                //printf("\tresults[%lu] = %lu\n", i, results[i].first);

                // Take unique_ptr
                term_ptr.swap(results[i].second);
                results.erase(results.begin() + i);
                break;
            }
        }

        if (!term_ptr.get()) {
            // Didn't find start; release lock and wait for worker_thread
            work_ready.wait(guard);
            continue;
        }

        // Release lock while doing testing
        assert(guard);
        guard.unlock();

        for (auto t : *term_ptr) {
            print_match(t.first, t.second);
        }

        // Overlap by segments by MAX_DIST, because easier than saving final X values.
        start += SEGMENT - MAX_DIST;

        // Relock so that loop starts with lock
        guard.lock();
    }

    // Do something non-standard for last segment so we end EXACTLY [start, stop)
    //      want to check each number, num, in range [start, stop) to see if num + d is a pair.
    print_results(start-1);
}


int main(int argc, char** argv) {
    // Allow comma seperators
    setlocale(LC_NUMERIC, "");

    //printf("Compiled with GMP %d.%d.%d\n",
    //    __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL);

    // Range is [START, STOP)

    uint64_t SEGMENT = 360360; //+ 11; //1 << 17;
    uint64_t START = 0;
    uint64_t STOP = 1e9;

    if (argc == 2)
        STOP = atol(argv[1]);

    if (START > 0)
        // Round to next lower segment
        START = ((START-1) / SEGMENT + 1) * SEGMENT;

    printf("Checking [%lu, %lu) in sets of %lu\n\n",  START, STOP, SEGMENT);

    assert(sigmaSelfCheck());

    A217259 runner(START, STOP, SEGMENT);

    // For single-threaded
    runner.iterate();

    // For multi-threaded
    //runner.multithreaded_iterate();
}
