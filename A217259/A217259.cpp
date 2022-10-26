// g++ -g -O3 --std=c++17 -Werror -Wall A217259.cpp -lgmpxx -lgmp -pthread -fopenmp && time ./a.out

/*

Results From [0, 9999999999] in groups of 120120
Summation: 17683099069429043286 (low64 bits)
Found 455,052,511 primes.
Found prime pairs at these distances:
	2 -> 27,412,679 pairs
	4 -> 27,409,999 pairs
	6 -> 54,818,296 pairs
	8 -> 27,411,508 pairs
	10 -> 36,548,839 pairs
Found 31 total composites, from these distances:
    3,1,11,3,11,2


Results From [0, 99999999999] in groups of 360360	took 219 seconds
Summation: 15869221272870986088 (low64 bits)
Found 4,118,054,813 primes.
Found prime pairs at these distances:
	2 -> 224,376,048 pairs
	4 -> 224,373,161 pairs
	6 -> 448,725,003 pairs
	8 -> 224,365,334 pairs
	10 -> 299,140,330 pairs
Found 31 total composites, from these distances:
    3,1,11,3,11,2


Results From [0, 999999999999] in groups of 360360	took 3257 seconds
Summation: 502152443576965830 (low64 bits)
Found 37,607,912,018 primes.
Found prime pairs at these distances:
	2 -> 1,870,585,220 pairs
	4 -> 1,870,585,459 pairs
	6 -> 3,741,217,498 pairs
	8 -> 1,870,580,394 pairs
	10 -> 2,494,056,601 pairs
Found 33 total composites, from these distances:
    3,1,12,3,12,2


Old result on 11 threads with possibly buggy code

225569500 	100,580,280,191 		150.0 seconds elapsed 670.5M/s
434699100 	205,449,033,509 		340.0 seconds elapsed 604.3M/s
625725600 	305,209,486,511 		540.0 seconds elapsed 565.2M/s
994193100 	504,381,603,587 		980.0 seconds elapsed 514.7M/s
1875582200	1,002,891,288,197		2260.0 seconds elapsed 443.8M/s
3556689700	2,002,379,757,239		5400.0 seconds elapsed 370.8M/s
8316959500	5,002,892,027,909		18100.0 seconds elapsed 276.4M/s
15838597700	10,002,668,115,959		46320.0 seconds elapsed 215.9M/s
30198403100	19,999,673,308,007		122420.0 seconds elapsed 163.4M/s

Old result on 11 threads with MAX_DIST=11 (and considers odds)

244550000 	107,630,754,461 		220.0 seconds elapsed 489.2M/s (interval: 491.4M/s)
460490000 	204,200,984,681 		420.0 seconds elapsed 486.2M/s (interval: 475.3M/s)
1159505000	509,224,217,501 		1080.0 seconds elapsed 471.5M/s (interval: 464.3M/s)
2327115000	1,009,483,795,889		2200.0 seconds elapsed 458.9M/s (interval: 445.6M/s)
4688615000	2,008,524,931,949		4640.0 seconds elapsed 432.9M/s (interval: 406.5M/s)
11857720000	5,006,181,689,021		11960.0 seconds elapsed 418.6M/s (interval: 400.4M/s)
23905255000	10,005,373,904,609		24840.0 seconds elapsed 402.8M/s (interval: 375.1M/s)
36014615000	15,012,087,118,661		38520.0 seconds elapsed 389.7M/s (interval: 359.5M/s)
40000000000	16,657,742,635,487 		twin prime
48130080000	20,012,474,218,211		52800.0 seconds elapsed 379.0M/s (interval: 343.6M/s)
60233975000	25,002,824,353,217		67600.0 seconds elapsed 369.9M/s (interval: 329.7M/s)
61877085000	25,679,962,609,589		69680.0 seconds elapsed 368.5M/s (interval: 317.2M/s)


twin prime count (up to 2e13) confirmed by primesieve in 432 seconds.
*/

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <mutex>
#include <thread>
#include <tuple>
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
const vector<uint64_t> SegmentedSieveOfSigma(uint64_t start, uint64_t N) {
    auto sums = vector<uint64_t>(N, 0);

    /* // Adjust to include n & 1 as a divisor of n
    for (uint64_t i = 0; i < N; i++) {
        sums[i] += start + i + 1;
    }
    // */

    if (start == 0) {
        sums[0] = -3;  // undefined
        sums[1] = -5;  // different undefined
    }

    auto past = start + N;
    uint64_t isqrt = calc_isqrt(past - 1);
    assert( isqrt * isqrt < past );
    assert( (isqrt+1) * (isqrt+1) >= past );

    // Handle the new squared factors (start <= f^2 < start+N) seperately
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

    uint64_t start_m_1 = start == 0 ? 0 : (start - 1);
    // Keeps track of what factors/divisors we've summed up.
    uint64_t factor = 2;

    /**
     * For very small factors break sieve_length into even smaller cached chunk.
     * This givese even better memory perforance
     */
    if (0) {
        constexpr uint32_t CHUNK_SIZE = 1 << 15; // 32768
        assert(N > 2 * CHUNK_SIZE);

        // [interval_index, interval_end)
        uint32_t i_i = 0;
        uint32_t i_end = CHUNK_SIZE;
        for (; i_i < N; i_i = i_end, i_end += CHUNK_SIZE) {
            i_end = std::max<uint32_t>(i_end, N);

            const uint64_t start_m_1 = start + i_i - 1;

            // Could possible cache count[factor] to avoid 2-10 recomputations
            for (factor = 2; factor <= std::min(isqrt, 1000ul); factor++) {
                uint64_t count = start_m_1 / factor + 1;
                uint32_t index = count * factor - start;
                uint64_t add = factor + count;

                assert(i_i <= index && index < i_end);
                assert( (i_end == N) || (index + 16*factor < i_end));

                /**
                 * NOTE: Made many previous attempts at loop optimization
                 * FIXME: Try again not that has better caching
                 * with #praga GCC unroll X
                 * with manual unrolling
                 * with 16x, 8x, 4x
                 * Have failed to make this sigificantly faster
                 */

                // Loop unrolled 4x for small factors
                for (; index + (factor<<2) < i_end; index += (factor<<2), add += 4) {
                    sums[index           ]  += add;
                    sums[index +   factor]  += add+1;
                    sums[index + 2*factor]  += add+2;
                    sums[index + 3*factor]  += add+3;
                }

                for (; index < N; index += factor, add++)
                    sums[index] += add;
            }
        }
    }


    for (; factor <= std::min(isqrt, N/5); factor++) {
        uint64_t count = start_m_1 / factor + 1;
        uint32_t index = count * factor - start;
        uint64_t add = factor + count;

        // Loop unrolled 4x for small factors
        for (; index + (factor<<2) < N; index += (factor<<2), add += 4) {
            sums[index           ]  += add;
            sums[index +   factor]  += add+1;
            sums[index + 2*factor]  += add+2;
            sums[index + 3*factor]  += add+3;
        }

        for (; index < N; index += factor, add++)
            sums[index] += add;
    }

    // Handles factors that can appear more than once
    for (; factor <= std::min(isqrt, N); factor++) {
        uint64_t count = start_m_1 / factor + 1;
        uint32_t index = count * factor - start;
        uint64_t add = factor + count;

        for (; index < N; index += factor, add++)
            sums[index] += add;
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

vector<uint32_t> gen_primes(uint32_t min, uint32_t n) {
    vector<bool> status(n, 1);

    vector<uint32_t> primes = {};
    if (2 >= min) {
        primes.push_back(2);
    }

    uint32_t p = 3;
    for (; p*p < n; p += 2) {
        if (status[p]) {
            if (p >= min)
                primes.push_back(p);
            for (uint32_t m = p*p; m < n; m += p) {
                status[m] = 0;
            }
        }
    }
    for (; p < n; p += 2) {
        if (status[p]) {
            if (p >= min)
                primes.push_back(p);
        }
    }
    return primes;
}

// TODO consider cache libdivide?
uint32_t count_factor(uint64_t n, uint32_t f) {
    uint32_t count = 0;
    while (n) {
        uint64_t d = n / f;
        if (d * f != n)
            return count;

        count++;
        n = d;
    }
    return count;
}

class SegmentedPrimeSieveSigma {
    public:
        SegmentedPrimeSieveSigma(uint64_t SEGMENT)
            : sieve_length(SEGMENT) {

            base_factored.resize(sieve_length);
            std::fill_n(&base_factored[0], sieve_length, std::make_pair(1ul,1ul));

            for (auto p : _primes) {
                if (sieve_length % p != 0)
                    continue;

                uint32_t max_power = count_factor(sieve_length, p);

                vector<uint64_t> powers = {1};
                for (auto power = uint32_t(1); power <= (max_power + 1); power++) {
                    powers.push_back(p * powers.back());
                }

                auto gcd = powers[max_power];
                assert(sieve_length % gcd == 0);

                sieve_primes.push_back(p);
                sieve_exponents.push_back(max_power);
                sieve_factors.push_back(gcd);

                // index=0 is WEIRD, handled manually in `next(...)`.

                // compute: 1 + p + p^2 + ...
                for (uint32_t i = p; i < sieve_length; i += p) {
                    // For all k what is the largest power of p to divide
                    // gcd(k * sieve_length + i, p^inf)

                    auto count = std::min(max_power, count_factor(i, p));
                    //base_factor[i] *= powers[count];
                    //base_product[i] *= (powers[count + 1] - 1) / (p - 1);
                    base_factored[i].first *= powers[count];
                    base_factored[i].second *= (powers[count + 1] - 1) / (p - 1);
                }
            }

            uint64_t p_mult = 1;
            for (auto g : sieve_factors)
                p_mult *= g;
            assert(p_mult == sieve_length);

            // Remove any entry in primes that's in sieve_primes
            {
                uint32_t i = 0;
                vector<uint32_t> filtered_primes;
                filtered_primes.reserve(_primes.size());
                for(auto p: _primes) {
                    if (sieve_primes[i] != p) {
                        filtered_primes.push_back(p);
                    } else {
                        if (i+1 < sieve_primes.size())
                            i++;
                    }
                }
                assert(filtered_primes.size() + sieve_primes.size() == _primes.size());
                _primes.swap(filtered_primes);
            }

            uint32_t cur = 0;
            for (auto p : _primes) {
                prime_gaps.push_back(p - cur);
                cur = p;
            }
        }

        const vector<uint64_t> next(uint64_t start);

        const uint32_t sieve_length;

    private:
        // More primes than we'll need so that list is constant
        vector<uint32_t> _primes = gen_primes(2, 7'100'000);
        // Stored as differences to reduce size
        vector<uint8_t> prime_gaps;

        // Generated from sieve_length
        vector<uint32_t> sieve_primes;
        vector<uint32_t> sieve_exponents;
        vector<uint32_t> sieve_factors;

        // Pre-handle sieve_factors
        vector<std::pair<uint64_t,uint64_t>> base_factored;
        //vector<uint64_t> base_factor;
        //vector<uint64_t> base_product;
};


const vector<uint64_t> SegmentedPrimeSieveSigma::next(uint64_t start) {
    assert(start % sieve_length == 0);
    assert((uint64_t) _primes.back() * _primes.back() > start + sieve_length);

    auto factored = vector<std::pair<uint64_t,uint64_t>>(base_factored);

    // build up (start + i) as p_1^e_1 * p_2^e_2 * ....
    //auto factors = vector<uint64_t>(base_factor);

    // Products of p_1^(e_1+1) / (p_1 - 1) over all p_i
    //auto products = vector<uint64_t>(base_product);


    uint32_t max_pp = 0;
    uint64_t prime_powers[32] = {1}; // 3^30 > 1e14, 30+1+1 = 32

    // Specialized p=2 loop
    {
        uint64_t p = 2;
        assert(p == sieve_primes[0]);
        uint32_t sieve_p = sieve_factors[0];
        uint32_t sieve_power = sieve_exponents[0];
        assert(sieve_length % p == 0);
        assert(sieve_length % sieve_p == 0);

        // V2, do multiples of p^6 (skipping multiples of p^7), p^7, p^8, ...
        uint64_t processing_power = sieve_power;
        uint64_t processing = sieve_p;
        uint64_t next_power = processing << 1;
        while (true) {
            processing = next_power;
            processing_power += 1;
            next_power <<= 1;

            uint64_t new_mult = (next_power - 1);

            // First multiple of processing AFTER start (to avoid weirdness as 0)
            uint64_t count = (start + processing) >> processing_power;
            uint64_t first_index = count * processing - start;
            if (first_index >= sieve_length)
                break;

            if ((count & 1) == 0) {
                // Multiple of next_power, skip
                first_index += processing;
                count += 1;
            }

            // Only need to fix indexes where count % p == 0 -> p, 2p, 3p, ...
            for (uint64_t index = first_index;
                    index < sieve_length;
                    count += 2, index += (processing<<1)) {

                //assert((count & 1) == 1);
                auto old_power = std::min(count_factor(index, p), sieve_power);
                factored[index].first <<= processing_power - old_power;

                uint64_t old_mult = (1ul << (old_power + 1)) - 1;
                factored[index].second = (factored[index].second / old_mult) * new_mult;

            }
        }

        // base[0] is WEIRD
        {
            uint64_t count = start / sieve_p;
            auto new_power = count_factor(count, p) + sieve_power;
            assert(factored[0].first % p != 0);
            factored[0].first <<= new_power;

            uint64_t new_mult = (1ul << (new_power + 1)) - 1;
            factored[0].second *= new_mult;
        }
    }

    for(uint32_t p_i = 1; p_i < sieve_primes.size(); p_i++) {
        uint64_t p = sieve_primes[p_i];
        uint32_t sieve_power = sieve_exponents[p_i];
        uint32_t sieve_p = sieve_factors[p_i];

        assert(sieve_length % sieve_p == 0);
        assert(sieve_length % p == 0);

        {
            for (uint32_t exp = 1; exp <= sieve_power; exp++) {
                prime_powers[exp] = p * prime_powers[exp - 1];
            }
            assert(prime_powers[sieve_power] == sieve_p);

            // need one extra for correcting products[i] with existing = sieve_power
            prime_powers[sieve_power+1] = p * prime_powers[sieve_power];
            max_pp = sieve_power+1;
        }

        // All numbers with fewer than sieve_power count of p are handled by base!

        /**
         * V1 Handled all multiples of p^(sieve_power+1) in a row
         * This required calculating `count_factor(index, p)`
         * and doesn't give better caching (remember p is large)
         * so that code has been removed.
         */

        // V2, do multiples of p^6 (skipping multiples of p^7) then multiples of p^7, ...
        uint32_t processing_power = sieve_power;
        uint64_t processing = sieve_p;
        uint64_t next_power = processing * p;
        while (true) {
            processing = next_power;
            processing_power += 1;
            next_power *= p;

            uint64_t new_mult = (next_power - 1) / (p-1);

            prime_powers[processing_power + 1] = next_power;
            max_pp += 1;

            // First multiple of processing AFTER start (to avoid weirdness as 0)
            uint64_t count = (start + processing) / processing;
            uint64_t first_index = count * processing - start;
            if (first_index >= sieve_length)
                break;

            uint64_t mod_p = count % p;

            // Only need to fix indexes where count % p == 0 -> p, 2p, 3p, ...
            for (uint64_t index = first_index;
                    index < sieve_length;
                    count += 1, index += processing) {
                //assert( (start + index) % processing == 0 );

                if (mod_p != 0) {
                    // new_power = processing_power
                    auto old_power = std::min(count_factor(index, p), sieve_power);
                    //assert(1 <= old_power);
                    //assert(old_power <= sieve_power);

                    //assert(factored[index].first % prime_powers[old_power] == 0);
                    factored[index].first *= prime_powers[processing_power - old_power];

                    uint64_t old_mult = (prime_powers[old_power + 1] - 1) / (p-1);
                    factored[index].second = (factored[index].second / old_mult) * new_mult;
                }

                mod_p += 1;
                if (mod_p == p)
                    mod_p -= p;
            }
        }

        // base[0] is WEIRD
        {
            uint64_t count = start / sieve_p;
            auto new_power = count_factor(count, p) + sieve_power;

            while (new_power >= max_pp) {
                prime_powers[max_pp + 1] = p * prime_powers[max_pp];
                max_pp += 1;
            }

            assert(max_pp > new_power);
            // No already added factors of p to "remove"
            uint64_t new_mult = (prime_powers[new_power + 1] - 1) / (p-1);
            factored[0].first *= prime_powers[new_power];
            factored[0].second *= new_mult;
        }
    }

    auto past = start + sieve_length;
    auto start_m_1 = start == 0 ? (2 - 1) : (start - 1);

    uint32_t p_i = 0;
    uint64_t p = 0;
    for (p_i = 0; ; p_i++) {
        p += prime_gaps[p_i];
        assert (p > 5);

        // Only one multiple in segment, handle without a loop below
        if (p >= sieve_length)
            break;

        // Only need to handle up to sqrt
        const uint64_t p2 = p*p;
        if (p2 >= past)
            break;

        // First multiple of p >= start
        uint64_t count = (start_m_1 + p) / p;
        uint32_t index = p * count - start;
        if (index >= sieve_length)
            continue;

        prime_powers[1] = p;
        prime_powers[2] = p2;
        prime_powers[3] = p2*p;
        prime_powers[4] = p2*p2;
        max_pp = 4;

        // Break loop over index into two parts
        //   indexes with factor of a single p
        //   indexes divisible by p^2 (or greater)
        // Do p - mod_pp iterations to get to next multiple of p^2

        uint64_t mod_pp = count % p;
        if (mod_pp != 0) {
            uint64_t simple_loops = (p - mod_pp);
            assert(mod_pp + simple_loops == p);

            // Update these for after the loop
            mod_pp = 0;
            count += simple_loops;

            for (; index < sieve_length && simple_loops > 0; index += p, simple_loops--) {
                // only a single power of p
                factored[index].first *= p;
                factored[index].second *= (1 + p);
            }
        }

        assert(mod_pp == 0);
        while (index < sieve_length) {
            { // 2 or more powers of pp (VERY RARE)
                auto pp = count_factor(count, p) + 1;
                while (pp >= max_pp) {
                    prime_powers[max_pp+1] = p * prime_powers[max_pp];
                    max_pp += 1;
                }
                factored[index].first *= prime_powers[pp];
                factored[index].second *= (prime_powers[pp] * p - 1) / (p - 1);

                index += p;
                count += 1;
            }

            // Do up to p-1 iterations of easy code
            count += p-1;

            // Maybe loop unroll 2x?
            // smaller of sieve_index and index_p[p-1]
            uint32_t stop_index = std::min<uint64_t>(sieve_length, index + p * (p-1));
            for (; index < stop_index; index += p) {
                factored[index].first *= p;
                factored[index].second *= (1 + p);
            }
        }
    }

    assert(p_i > 0);
    for (; ; p_i++) {
        p += prime_gaps[p_i];
        const uint64_t p2 = p*p;
        if (p2 >= past)
            break;

        // First multiple of p >= start
        uint64_t count = (start_m_1 + p) / p;
        uint32_t index = p * count - start;

        if (index < sieve_length) {
            // Check if count is a multiple of p

            if (count != (count / p) * p) {
                [[likely]]
                // Common (and easy) case
                factored[index].first *= p;
                factored[index].second *= (1 + p);
            } else {
                // Hard case - More than one factor of p - VERY RARE
                auto pp = count_factor(count, p) + 1;
                uint64_t prime_power = p2;
                for(; pp > 2; pp--) {
                    prime_power *= p;
                }
                factored[index].first *= prime_power;
                factored[index].second *= (prime_power * p - 1) / (p - 1);
            }
        }
    }

    uint64_t num = start;

    vector<uint64_t> results(sieve_length);

    assert(sieve_length % 6 == 0);
    #pragma GCC unroll 6
    #pragma GCC ivdep
    for (uint32_t i = 0; i < sieve_length; i ++, num ++) {
        uint64_t prime = num / factored[i].first;
        //assert(num == 0 || prime == 1 || ((prime*prime >= num) && ((prime&1) == 1)));
        //assert(num == 0 || num % prime == 0);
        results[i] = factored[i].second * ((prime > 1) + prime) - (num+1);
    }

    return results;
}

uint32_t calc_bucket_capacity(uint64_t stop, uint64_t min_f, uint64_t sieve_length) {
    // Sadly this is not my oldest friend Mertens, but is instead harmonic numbers
    // FIXME: Doesn't seem to be working in all cases? Maybe start / stop is affecting this?
    float max_divisor = std::sqrt(stop + 2 * sieve_length);
    float estimate_per_num = std::max(0.0, std::log(max_divisor) - std::log(min_f));
    uint32_t estimate = estimate_per_num * sieve_length * 1.05 + 10;
    printf("\tWith stop=%lu (%.0f) guessing bucket_capacity=%u\n", stop, max_divisor, estimate);
    printf("\t\tbuckets use starts at %'lu\n", min_f * min_f);
    return estimate;
}

class SegmentedBucketedSigma {
    public:
        SegmentedBucketedSigma(uint64_t start, uint64_t stop, uint64_t length) :
                sieve_length(length),
                MIN_BUCKET_DIVISOR(3*length),
                bucket_capacity(calc_bucket_capacity(stop, MIN_BUCKET_DIVISOR, sieve_length)) {
            //assert(sieve_length == 2*3*2*5*1*7*2*3*1*11*1*1*13*1*1*2);
            //assert(sieve_length == 360360);
            assert(sieve_length % (2*3*2*5) == 0);
            assert(100000 < sieve_length);
            assert(sieve_length < 1'800'000);

            base_counts.resize(sieve_length);
            base_delta.resize(sieve_length);

            // So that we don't have to deal with f^2 entries after 0th interval.
            assert(MAX_BASE_F*MAX_BASE_F < sieve_length);

            assert((1 << 5) == MAX_BUCKET + 1); // Must be 1 off a power of two
            assert(bucket_capacity < 4'000'000);
            for (auto &bucket : bucketed_divisors) {
                bucket.resize(bucket_capacity);
            }

            sums.resize(sieve_length);

            // Handles setting base_count, base_delta, buckets
            jump_to(start);
        }

        void jump_to(uint64_t new_start);
        const vector<uint64_t>& next(uint64_t verify_start);

    private:
        const uint64_t sieve_length;

        // Results are stored here!
        vector<uint64_t> sums;

        // Start of the next interval, must be a multiple of sieve_length
        uint64_t start = 0;
        // Largest factor added to bucket_divisors
        uint64_t next_factor = 0;

        // Don't bucket divisor less than X
        const uint64_t MIN_BUCKET_DIVISOR;

        // This limits STOP to ~ (MAX_BUCKET * SIEVE_LENGTH)^2 = 132e12
        // Needs to be a power of 2 minus 1
        static const uint32_t MAX_BUCKET = 31;
        // divisors above MIN_BUCKET_DIVISOR goes in bucket[(next_multiple / sieve_length) % (MAX_BUCKET+1)]
        vector<uint32_t> bucketed_divisors[MAX_BUCKET + 1];
        const uint64_t bucket_capacity;

        const uint32_t MAX_BASE_F = 300; //600; // isqrt(360360-1)

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
};


void SegmentedBucketedSigma::jump_to(uint64_t new_start) {
    assert(new_start % sieve_length == 0);
    start = new_start;

    uint64_t current_bucket = start / sieve_length;

    std::fill(base_counts.begin(), base_counts.end(), 0);
    std::fill(base_delta.begin(), base_delta.end(), 0);

    // Number of factors handled "automatically"
    float handled_by_base = 0;
    for (uint32_t f = 2; f <= MAX_BASE_F; f++) {
        if (sieve_length % f == 0) {
            handled_by_base += (sieve_length / f);

            // Account for base_delta being added once
            for (uint32_t i = 0; i < sieve_length; i += f)
                base_counts[i] += f + (start + i) / f;

            for (uint32_t i = 0; i < sieve_length; i += f)
                base_delta[i] += (sieve_length / f);
        }
    }
    printf("\tSet up sieve of length %lu handles %.0f automatically = %.2f/num\n",
            sieve_length, handled_by_base, handled_by_base / sieve_length);

    // Spot check
    assert(base_counts[2] == 2 + (new_start + 2) / 2);
    assert(base_counts[3] == 3 + (new_start + 3) / 3);
    assert(base_counts[5] == 5 + (new_start + 5) / 5);

    for (auto &bucket : bucketed_divisors) {
        bucket.clear();
    }

    // Everything [MIN_BUCKET_DIVISOR, isqrt(start)) will be added to a bucket.
    next_factor = start == 0 ? 0 : calc_isqrt(start - 1) + 1;

    for (uint64_t f = MIN_BUCKET_DIVISOR; f < next_factor; f++) {
        assert(f * f < start);
        /**
         * Distance, from start, to next multiple of f
         * This is conceptually -start % f
         */
        uint64_t count = (start + f - 1) / f;
        assert(count > f);
        uint64_t index = count * f;
        assert(index >= start);

        uint64_t bucket = index / sieve_length;
        assert(bucket >= current_bucket);
        assert((bucket-current_bucket) <= MAX_BUCKET);

        // & MAX_BUCKET is the same as % (MAX_BUCKET+1)
        bucketed_divisors[bucket & MAX_BUCKET].push_back(f);
    }
 }

const vector<uint64_t>& SegmentedBucketedSigma::next(uint64_t verify_start) {
    assert(start == verify_start);

    // TODO figure out how to relax this.
    // First interval is WEIRD with many f^2 < sieve_length
    if (start == 0) {
        auto t = SegmentedSieveOfSigma(0, sieve_length);
        for (uint64_t i = 0; i < sieve_length; i++) sums[i] = t[i];

        // let jump_to() correctly update state
        jump_to(start + sieve_length);
        return sums;
    }

    for (uint64_t i = 0; i < sieve_length; i++) {
        sums[i] = base_counts[i];
        base_counts[i] += base_delta[i];
    }

    uint64_t current_bucket = start / sieve_length;

    // Hand unrolled loops for very small factor
    uint64_t factor = 2;
    for (; factor < std::min<uint32_t>(sieve_length / 9, next_factor); factor++) {
        if (factor <= MAX_BASE_F && sieve_length % factor == 0)
            continue;

        uint64_t count = (start + factor - 1) / factor;
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
        for (; updates >= 8; updates -= 8, add += 8, index += factor<<3) {
            sums[index           ]  += add;
            sums[index +   factor]  += add+1;
            sums[index + 2*factor]  += add+2;
            sums[index + 3*factor]  += add+3;
            sums[index + 4*factor]  += add+4;
            sums[index + 5*factor]  += add+5;
            sums[index + 6*factor]  += add+6;
            sums[index + 7*factor]  += add+7;
        }

        for (; updates > 0; updates--, add++, index += factor) {
            sums[index] += add;
        }

        assert((uint32_t) index >= sieve_length);
    }

    assert(MIN_BUCKET_DIVISOR > sieve_length);

    // Handles factors that appear at least once (but possible more times)
    for (; factor < std::min(sieve_length, next_factor); factor++) {
        uint64_t count = (start + factor - 1) / factor;
        uint32_t index = count * factor - start;
        auto add = count + factor;

        for (; index < sieve_length; index += factor, add++)
            sums[index] += add;
    }

    // Handles larger factors that can only appear once but aren't bucketed
    for (; factor < std::min(next_factor, MIN_BUCKET_DIVISOR); factor++) {
        uint64_t count = (start + factor - 1) / factor;
        uint32_t index = count * factor - start;
        if (index < sieve_length) {
            sums[index] += count + factor;
        }
    }

    // Handles larger factors in this bucket.
    for (uint64_t factor : bucketed_divisors[current_bucket & MAX_BUCKET]) {
        assert(factor >= MIN_BUCKET_DIVISOR);
        uint64_t count = (start + factor - 1) / factor;
        uint64_t mult = count * factor;
        //assert(mult >= start);
        uint64_t index = mult - start;
        //assert(index < sieve_length);
        sums[index] += count + factor;

        // TODO libdivide for "/ sieve_length"
        // add to a new bucket
        mult += factor;
        uint64_t bucket = mult / sieve_length;
        //assert((bucket-current_bucket) <= MAX_BUCKET);
        bucketed_divisors[bucket & MAX_BUCKET].push_back(factor);
    }
    {
        auto size = bucketed_divisors[current_bucket & MAX_BUCKET].size();
        if (size > bucket_capacity || (size > 0 && start % 1001 == 1))
            printf("%lu pulled %lu from current_bucket! sieve_length: %lu, bucket_capacity: %lu \n",
                    start, size, sieve_length, bucket_capacity);
        if (size > bucket_capacity)
            exit(1);
    }

    // Clear all the processed factors in buckets.
    bucketed_divisors[current_bucket & MAX_BUCKET].clear();

    // New squared factors
    {
        // Otherwise have to not count these.
        assert(next_factor > MAX_BASE_F);

        auto past = start + sieve_length;
        uint64_t f = next_factor;
        for (; ; f++) {

            uint64_t f_2 = f * f;
            assert(f_2 >= start);

            if (f_2 >= past)
                break;

            uint64_t index = f_2 - start;
            assert(index < sieve_length);

            // Update the square in this segmant/interval.
            sums[index] += f;

            if (f >= MIN_BUCKET_DIVISOR) {
                // Add next multiple to a bucket
                uint64_t bucket = (f_2 + f) / sieve_length;
                assert((bucket-current_bucket) <= MAX_BUCKET);
                bucketed_divisors[bucket & MAX_BUCKET].push_back(f);
            } else {
                // Maybe have to update more multiples in this segmant.
                index += f;
                auto add = f + (f + 1);
                for (; index < sieve_length; index += f, add ++) {
                    sums[index] += add;
                }
            }
        }
        next_factor = f;
    }

    start += sieve_length;
    return sums;
}

class SegmentedSieveSigma {
    public:
        SegmentedSieveSigma(uint64_t start, uint64_t length) : sieve_length(length) {
            //assert(sieve_length == 2*3*2*5*1*7*2*3*1*11*1*13*1*1);
            //assert(sieve_length == 360360);

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

        const uint32_t MAX_BASE_F = 407; // isqrt(360360-1)

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
    assert(base_counts[2] == 2 + (new_start + 2) / 2);
    assert(base_counts[3] == 3 + (new_start + 3) / 3);
    assert(base_counts[5] == 5 + (new_start + 5) / 5);
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
        //assert(count != -1ul);
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
        for (; updates >= 8; updates -= 8, add += 8, index += factor<<3) {
            sums[index           ]  += add;
            sums[index +   factor]  += add+1;
            sums[index + 2*factor]  += add+2;
            sums[index + 3*factor]  += add+3;
            sums[index + 4*factor]  += add+4;
            sums[index + 5*factor]  += add+5;
            sums[index + 6*factor]  += add+6;
            sums[index + 7*factor]  += add+7;
        }

        for (; updates > 0; updates--, add++, index += factor) {
            sums[index] += add;
        }

        assert((uint32_t) index >= sieve_length);
        counts[factor] = add - factor;
    }

    // Handles factors that appear at least once (but possible more times)
    for (; factor < std::min<uint32_t>(counts.size(), sieve_length); factor++) {
        auto count = counts[factor];
        //assert(count != -1ul);
        //assert(count  == std::max(factor + 1, ((start + factor - 1) / factor)));
        uint32_t index = count * factor - start;
        auto add = count + factor;

        for (; index < sieve_length; index += factor, add++)
            sums[index] += add;

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
    auto errors = 0;

    // Several tests
    for (const auto [START, N, EXPECTED] : {
            std::tuple{0ul, 4'000'000ul, 13159467448256ul},
            {(((uint64_t) 5e12) / 1441440) * 1441440, 1'000'000, 3224658678984463479ul},
            {(((uint64_t) 12e12) / 1441440) * 1441440, 100'000, 773889434291175054}}) {

        const uint64_t STOP = START + N;

        printf("\n\tSelf Check [%lu, %lu) sum: %lu\n", START, STOP, EXPECTED);

        vector<uint64_t> sigmas_a, sigmas_b, sigmas_c, sigmas_d;
        {
            for (uint64_t start = START; start < STOP; start += 100'000) {
                auto temp = SegmentedSieveOfSigma(start, 100'000);
                sigmas_a.insert(std::end(sigmas_a), std::begin(temp), std::end(temp));
            }
        }
        {
            const uint64_t SEGMENT = 360360;
            auto sieveSum  = SegmentedSieveSigma(START, SEGMENT);
            for (uint64_t start = START; start < STOP; start += SEGMENT) {
                auto temp = sieveSum.next(start);
                sigmas_b.insert(std::end(sigmas_b), std::begin(temp), std::end(temp));
            }
        }
        {
            const uint64_t SEGMENT = 360360;
            auto sieveProd = SegmentedPrimeSieveSigma(SEGMENT);
            for (uint64_t start = START; start < STOP; start += SEGMENT) {
                auto temp = sieveProd.next(start);
                sigmas_c.insert(std::end(sigmas_c), std::begin(temp), std::end(temp));
            }
        }
        {
            const uint64_t SEGMENT = 110880;
            auto sieveBucket = SegmentedBucketedSigma(START, STOP+SEGMENT, SEGMENT);
            for (uint64_t start = START; start < STOP; start += SEGMENT) {
                auto temp = sieveBucket.next(start);
                sigmas_d.insert(std::end(sigmas_d), std::begin(temp), std::end(temp));
            }
        }

        assert(sigmas_a.size() >= N);
        assert(sigmas_b.size() >= N);
        assert(sigmas_c.size() >= N);
        assert(sigmas_d.size() >= N);

        uint64_t summation = 0;
        for (uint32_t i = 2; i < N; i++) {
            summation += sigmas_a[i] + i + 1;

            if (sigmas_a[i] != sigmas_b[i] ||
                    sigmas_a[i] != sigmas_c[i] ||
                    sigmas_a[i] != sigmas_d[i] ) {
                printf("Self Check mismatch at %lu (%u) | %lu, %lu, %lu, %lu\n",
                        START + i, i,
                        sigmas_a[i], sigmas_b[i],
                        sigmas_c[i], sigmas_d[i]);
                errors += 1;
                if (errors > 5)
                    return false;
            }
        }

        if (summation != EXPECTED) {
            printf("\tSummation %lu != %lu Self-Check Failed\n", summation, EXPECTED);
            assert(summation == EXPECTED);
        }
    }
    printf("\n");
    return errors == 0;
}

class A217259 {
    public:
        A217259(uint64_t start, uint64_t stop, uint64_t segment_size)
            : START(start), STOP(stop), SEGMENT(segment_size),
              SEGMENTED_START((START / SEGMENT) * SEGMENT)
        {
            // Need to add 1 + n for all n in range(2, STOP)
            sum_sigma = (STOP*STOP + STOP - 6) / 2;

            print_interval_i = SEGMENTED_START;
        }


        bool check_known_composite_match(int16_t dist, uint64_t i);
        void print_match(int16_t dist, uint64_t i);
        bool print_match_and_check_known(int16_t dist, uint64_t i) {
            print_match(dist, i);
            if (dist < 0)
                // All prime gaps are "known"
                return true;
            return check_known_composite_match(dist, i);
        }
        void print_results(uint64_t last_included, uint64_t summation);

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
        // First multiple of SEGMENT that contains START
        const uint64_t SEGMENTED_START;

        static const uint32_t MAX_DIST = 10;

        uint64_t sum_sigma = 0;

        int64_t found_prime[MAX_DIST+1] = {};
        int64_t total_composites = 0;
        int64_t found_composite[MAX_DIST+1] = {};
        int64_t print_mult = 1;
        uint64_t next_time = 5;
        uint64_t print_interval_i = 0;

        // time for START and INTERVAL
        std::chrono::time_point<std::chrono::system_clock> t_START =
            std::chrono::system_clock::now();
        std::chrono::time_point<std::chrono::system_clock> t_INTERVAL =
            t_START;

        // Guard for results
        std::mutex g_control;
        std::condition_variable work_ready;

        // Multithreaded results output
        typedef vector<std::pair<int8_t,uint64_t>> terms_t;
        // START, start&end of segment, matches inside segment
        vector<std::tuple<uint64_t,vector<uint64_t>*,terms_t*>> results;
};


void A217259::print_results(uint64_t last_stop, uint64_t summation) {
    double elapsed = std::chrono::duration<double>(std::chrono::system_clock::now() - t_START).count();

    printf("\n\n\n\n");
    printf("Results From [%'lu,  %'lu] in groups of %lu\ttook %.1f seconds",
            START, last_stop, SEGMENT, elapsed);

    // TODO track if we ran singlethreaded or multithreaded
    printf("\n");
    if (summation)
        printf("Summation: %lu (low64 bits)\n", summation);
    printf("Found %'lu primes.\n", found_prime[0]);
    printf("Found prime pairs at these distances:\n");
    for (uint32_t d = 1; d <= MAX_DIST; d++) {
        if (found_prime[d] > 1) // 2 and N isn't iteresting
            printf("\t%d -> %'lu pairs\n", d, found_prime[d]);
    }
    printf("Found %lu total composites, from these distances:\n", total_composites);
    for (uint32_t d = 1; d <= MAX_DIST; d++) {
        if (found_composite[d])
            printf("\t%d -> %'lu pairs\n", d, found_composite[d]);
    }
}

bool A217259::check_known_composite_match(int16_t dist, uint64_t i) {
    /**
     * Already known terms
     * DIST = 2 (A217259/A050507), 6 (A054903), 7 (A063680), 8 (A059118), 12 (A054902)
     * Smallest even terms: A054905
     * For odd terms see README.md and check_squares.py
     */
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

    if (dist == 13)
        for (uint64_t known : std::initializer_list<uint64_t>{
                4418
                })
            if (i == known)
                return true;


    // use sgn(d) to check if this matches sigma result
    for (uint64_t test : {i, i + dist}) {
        mpz_class t = test;
        int primality = mpz_probab_prime_p(t.get_mpz_t(), 20);
        if (primality != 0) {
            printf("\n\nPRIMALITY MISMATCH FOR %lu | %d vs %d\n\n", test, primality, dist);
            exit(1);
        }

        if (primality == 0 && test == i) {
            printf("\n\n\n\nNEW! %lu, %lu+%u | %lu is composite(%d)!\n\n\n\n\n",
                    i, i, abs(dist), i, primality);
        }
    }
    return false;
}

void A217259::print_match(int16_t dist, uint64_t i) {
    // A composite!
    if (dist > 0) {
        total_composites += 1;
        found_composite[dist] += 1;
    } else {
        found_prime[-dist] += 1;
    }

    if (dist > 0) {
        printf("%2ld composite\t%'-16lu \t\t%ld'th composite with d=%d\n",
                total_composites, i, found_composite[dist], dist);
    } else if (dist == -2) {
        auto twin_primes = found_prime[2];

        // So we can verify against A146214
        if (twin_primes == 5*print_mult)
            print_mult *= 10;

        if (twin_primes % print_mult == 0) {
            printf("%-10ld\t%'-16lu \t\ttwin prime\n", twin_primes, i);
        } else if (twin_primes % 5000 == 0) {
            // Avoid calling sys_clock on every find.
            auto now = std::chrono::system_clock::now();
            double elapsed = std::chrono::duration<double>(now - t_START).count();

            if (elapsed > next_time) {
                double elapsed_interval = std::chrono::duration<double>(now - t_INTERVAL).count();
                float rate = (i - START) / elapsed / 1e6;
                float rate_interval = (i - print_interval_i) /elapsed_interval / 1e6;

                printf("%-10ld\t%'-16lu\t\t%.1f seconds elapsed %.1fM/s (interval: %.1fM/s)\n",
                        twin_primes, i, elapsed, rate, rate_interval);

                // 5,10,15 ... 95,100,110,120,130 ... 190,200,220,240 ... 960,1000,1080,1160
                next_time += 5 * (1 + (next_time >= 30)) * (1 + (next_time >= 200)) * (1 + (next_time >= 1000));
                print_interval_i = i;
                t_INTERVAL = now;
            }
        }
    }
}

__attribute__((always_inline))
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
            if (!print_match_and_check_known(d, a))
                printf("\n\tsigmas[%lu] = %lu, sigmas[%lu] = %lu\n",
                    a, sigma_a, a + d, sigma_b);
        }
    }
}


void A217259::iterate() {
    // EXCLUDES STOP

    uint64_t last_sigmas[MAX_DIST];
    for (size_t i = 0; i < MAX_DIST; i++) {
        last_sigmas[i] = -200ul;
    }

    //auto sieve = SegmentedSieveSigma(SEGMENTED_START, SEGMENT);
    //auto sieve = SegmentedBucketedSigma(SEGMENTED_START, STOP+SEGMENT, SEGMENT);
    auto sieve = SegmentedPrimeSieveSigma(SEGMENT);

    for (uint64_t start = SEGMENTED_START; start < (STOP + MAX_DIST); start += SEGMENT) {
        //auto sigmas = SegmentedSieveOfSigma(start, SEGMENT);
        auto sigmas = sieve.next(start);

        // for (uint32_t i = 0; i < SEGMENT; i++)
        //     if (sigmas[i] != sigmas2[i])
        //         printf("\n\nSIGMA MISMATCH @ %lu (%u) | %lu vs %lu\n\n\n",
        //                 start + i, i, sigmas[i], sigmas2[i]);

        // Considering keep sum of sigma as a spot check

        // test i less than or equal to these values in the three loops
        uint64_t last_i_head = MAX_DIST-1;
        uint64_t last_i_main = SEGMENT-MAX_DIST-1;
        uint64_t last_i_tail = SEGMENT-1;

        // Start from std::max(2, START, start);
        uint64_t i_start = std::max(std::max(2ul, START), start) - start;

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
            for (uint32_t i = 0; i <= last_i_head; i ++) {
                // d >= MAX_DIST - i
                // round d up to the next even number
                for (uint32_t d = MAX_DIST - (i - (i&1)); d <= MAX_DIST; d += 2) {
                    //assert(i + d >= MAX_DIST);
                    //assert((d & 1) == 0);
                    process_pair(d, start - MAX_DIST + i, last_sigmas[i], sigmas[i + d - MAX_DIST]);
                }
            }
        }

        // MAIN checking
        for (uint32_t i = i_start; i <= last_i_main; i++) {
            sum_sigma += sigmas[i];
            found_prime[0] += (sigmas[i] == 0);

            auto sigmas_i = sigmas[i];
            // +10-20% speed up from unrolling this loop
            // This loop represents a lot of the iterations of the program
            #pragma GCC unroll (MAX_DIST/2 - 1)
            for (uint32_t d = 2; d <= MAX_DIST; d += 2) {
                process_pair(d, start + i, sigmas_i, sigmas[i + d]);
            }
        }

        // Checking of final items which may roll over to next window.
        for (uint32_t i = SEGMENT - MAX_DIST; i <= last_i_tail; i++) {
            last_sigmas[i - (SEGMENT-MAX_DIST)] = sigmas[i];

            sum_sigma += sigmas[i];
            found_prime[0] += (sigmas[i] == 0);

            // All the pairs which stay inside [start, start + SEGMENT)
            for (uint32_t d = 2; d < (SEGMENT - i); d += 2) {
                process_pair(d, start + i, sigmas[i], sigmas[i + d]);
            }
        }
    }

    print_results(STOP-1, sum_sigma);
}


void A217259::worker_thread() {
    // Work horse of multithreaded_iterate() {

    // Broken because of use of state
    //auto sieve = SegmentedSieveSigma(SEGMENTED_START, SEGMENT);
    //auto sieve = SegmentedBucketedSigma(SEGMENTED_START, STOP+SEGMENT, SEGMENT);

    // Work in multithreaded context
    auto sieve = SegmentedPrimeSieveSigma(SEGMENT);

    // TODO tune dynamic 4
    #pragma omp parallel for schedule(dynamic, 4)
    for (uint64_t start = SEGMENTED_START; start < (STOP + MAX_DIST); start += SEGMENT) {
        // Calculate results
        //auto sigmas = SegmentedSieveOfSigma(start, SEGMENT);
        auto sigmas = sieve.next(start);

        // Create vectors that will be added to results
        terms_t *terms = new terms_t();
        auto *partial_sigma = new vector<uint64_t>(2*MAX_DIST);
        {
            // Place initial and final terms with special dist marker
            // This allows worker_coordinator to find matches over segment boundaries
            uint32_t j = 0;
            for (uint32_t i = 0; i < MAX_DIST; i++)
                (*partial_sigma)[j++] = sigmas[i];

            for (uint32_t i = (SEGMENT-MAX_DIST); i < SEGMENT; i++)
                (*partial_sigma)[j++] = sigmas[i];

            assert(partial_sigma->size() == 2*MAX_DIST);
        }


        // Start from std::max(2, START, start);
        uint64_t start_i = std::max(std::max(2ul, START), start) - start;

        auto stop_i = SEGMENT-MAX_DIST;
        if (start + stop_i > STOP) {
            if (start < STOP) {
                stop_i = STOP - start;
            } else {
                stop_i = 0;
            }
        }
        assert(stop_i + MAX_DIST <= SEGMENT);

        uint32_t found_primes = 0;
        uint64_t sum_sigma_local = 0;

        // Only look at i in the non-overlapping part
        for (uint32_t i = start_i; i < stop_i; i++) {
            sum_sigma_local += sigmas[i];
            found_primes += (sigmas[i] == 0);

            // +20-50% speed up from unrolling this loop!
            // This loop represents a lot of the iterations of the program
            /*
            #pragma GCC unroll (MAX_DIST/2 - 1)
            for (uint16_t dist = 2; dist <= MAX_DIST; dist += 2) {
                if (sigmas[i] == sigmas[i + dist]) {
                    if (sigmas[i] == 0 && sigmas[i+dist] == 0) {
                        // Prime pair (e.g. twin, cousin, ...)
                        terms->push_back({-dist, start + i});
                    } else {
                        if (!check_known_composite_match(dist, start + i))
                            printf("\n\tsigmas[%lu] = %lu, sigmas[%lu] = %lu\n",
                                start+i, sigmas[i], start+i+dist, sigmas[i+dist]);

                        terms->push_back({dist, start + i});
                    }
                }
            }
            */
            if (sigmas[i] == 0) {
                // Search adjacent for prime pair (e.g. twin, cousin, ...)
                #pragma GCC unroll (MAX_DIST/2 - 1)
                for (uint16_t dist = 2; dist <= MAX_DIST; dist += 2)
                    if (sigmas[i + dist] == 0)
                        terms->push_back({-dist, start + i});
            } else {
                [[likely]]
                #pragma GCC unroll (MAX_DIST/2 - 1)
                for (uint16_t dist = 2; dist <= MAX_DIST; dist += 2) {
                    if (sigmas[i] == sigmas[i + dist]) {
                        if (!check_known_composite_match(dist, start + i))
                            printf("\n\tsigmas[%lu] = %lu, sigmas[%lu] = %lu\n",
                                start+i, sigmas[i], start+i+dist, sigmas[i+dist]);

                        terms->push_back({dist, start + i});
                    }
                }
            }
        }

        // Wait till I can safely queue results.
        std::unique_lock<std::mutex> guard(g_control);
        assert(guard.owns_lock());

        sum_sigma += sum_sigma_local;
        found_prime[0] += found_primes;

        //printf("\t\tWork %lu -> %lu ready\n", start, terms.size());
        results.emplace_back(start, partial_sigma, terms);

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

    // Rollover terms from last segment
    vector<uint64_t> sigmas;

    uint64_t start = SEGMENTED_START;
    while (start < STOP + MAX_DIST) {
        if (!guard) {
            guard.lock();
        }

        if (results.empty()) {
            // Release lock and wait for worker_thread
            work_ready.wait(guard);
            continue;
        }

        vector<uint64_t> *partial_sigma_ptr = nullptr;
        terms_t *term_ptr = nullptr;

        // Check if start in results
        for (size_t i = 0; i < results.size(); i++) {
            if (std::get<0>(results[i]) == start) {
                //printf("\tresults[%lu] = %lu\n", i, results[i].first);

                // Get pointer to vector of partial sigma
                partial_sigma_ptr = std::get<1>(results[i]);
                // Get pointer to vector of results
                term_ptr = std::get<2>(results[i]);

                results.erase(results.begin() + i);
                break;
            }
        }

        if (term_ptr == nullptr) {
            // Didn't find start; release lock and wait for worker_thread
            work_ready.wait(guard);
            continue;
        }

        // Release lock while doing testing
        assert(guard);
        guard.unlock();

        assert(partial_sigma_ptr->size() == 2*MAX_DIST || start+SEGMENT >= STOP);

        uint64_t sum_sigma_local = 0;
        uint32_t found_primes = 0;
        if (start > START) {
            for (uint64_t i = 0; i < MAX_DIST; i++) {
                sigmas.push_back((*partial_sigma_ptr)[i]);
            }

            assert(sigmas.size() == 2*MAX_DIST);

            uint32_t last_i_head = MAX_DIST;
            if (STOP <= start) {
                last_i_head = STOP + MAX_DIST - start;
            }

            for (uint32_t i = 0; i < last_i_head; i++) {
                sum_sigma_local += sigmas[i];
                found_primes += (sigmas[i] == 0);

                // Bounds are hard, verify explicitly.
                uint64_t num = start - MAX_DIST + i;
                assert(START <= num && num < STOP);

                #pragma GCC unroll (MAX_DIST/2 - 1)
                for (uint32_t d = 2; d <= MAX_DIST; d += 2) {
                    // First from last segment (AKA i < MAX_DIST)
                    // Second from this segment (AKA i >= MAX_DIST)
                    /*
                    printf("%lu (%u) -> %lu | %lu (%u) -> %lu\n",
                            start - MAX_DIST + i, i, sigmas[i],
                            start - MAX_DIST + i + d, i + d, sigmas[i + d]);
                    // */
                    process_pair(d, num, sigmas[i], sigmas[i + d]);
                }
            }
        }

        // Save the last MAX_DIST terms for comparison with the next SEGMENT's first terms.
        sigmas.clear();
        for (uint64_t i = MAX_DIST; i < 2*MAX_DIST; i++) {
            // -(MAX_DIST-1), -(MAX_DIST-2), ... -1
            sigmas.push_back((*partial_sigma_ptr)[i]);
        }

        {
            guard.lock();
            sum_sigma += sum_sigma_local;
            found_prime[0] += found_primes;
            guard.unlock();
        }

        for (auto t : *term_ptr) {
            print_match(t.first, t.second);
        }

        delete partial_sigma_ptr;
        delete term_ptr;

        start += SEGMENT;
    }

    // Last segment is correctly handled in worker_thread()
    print_results(STOP-1, sum_sigma);
}


int main(int argc, char** argv) {
    // Allow comma seperators
    setlocale(LC_NUMERIC, "");

    //printf("Compiled with GMP %d.%d.%d\n",
    //    __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL);

    printf("Running Self Check (DISABLE IN VALGRIND!)\n");
    assert(sigmaSelfCheck());

    // Range is [START, STOP)
    uint64_t START = 20e12;
    uint64_t STOP  = 20.002e12;

    /**
     * Sieve length lets in many of the sieves "precompute" divisors that divide the length
     *
     * Quality is roughly Sum(1/d, d divides sieve_length) = sigma(sieve_length, -1)
     *
     * If this is too large threads step on each other (and more threads will harm cache perf)
     * 150000-250000 is maybe optimal on my Ryzen 3900x?
     * These values all have great cover (3.15+)
     * 55440, 110880, 166320, 277200, 332640, 360360, 554400, 665280, 720720
     * And good cover
     * 90720, 234460,
     *
     * Compare 166320 with 178080 to see -10% performance impact
     */
    //uint64_t SEGMENT = 360360;  // Great for SegmentedSieveOfSigma
    //uint64_t SEGMENT = 166320;
    //uint64_t SEGMENT = 720720;

    // For SegmentedPrimeSieveSigma try 60060, 120120, 180180, 240240 (best for me), 300300, 360360
    uint64_t SEGMENT = 4*60060;

    /**
     * On powersave core
     * | SIEVE                    | 5T    | 10T   | 20T   | 50T   |
     * | SegmentedSieveOfSigma    | 11.3  | 8.3   | 6.1   | 3.9   | <- roughly 2.5e7 / sqrt(STOP)
     * | SegmentedSieveSigma      | 29.9  | 25.3  | 21.0  | 15.8  | <- uses counts[divisor]
     * | SegmentedBucketedSigma   | 19.2  | 15.9  | 13.6  | 11.9   | <- using much larger SEGMENT
     *      SegmentedSieveOfSigma doesn't have buckets
     *      SegmentievedSieveSigma has counts
     * | SegmentedPrimeSieveSigma | 48.5  | 44.4  | 39.0  | 33.5  | <- using 240240
     */

    if (argc == 2) {
        START = 0;
        STOP = atol(argv[1]);
    }

    printf("Checking [%'lu,  %'lu) in sets of %lu\n\n",  START, STOP, SEGMENT);

    A217259 runner(START, STOP, SEGMENT);

    // For single-threaded
    //runner.iterate();

    // For multi-threaded
    runner.multithreaded_iterate();
}
