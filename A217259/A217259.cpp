// g++ -g -O3 --std=c++17 -Werror -Wall A217259.cpp -lgmpxx -lgmp -pthread -fopenmp && time ./a.out

/*

multi-threaded(4) 2022/10/13 results

103743400  43,072,515,048  		110.0 seconds elapsed 391.6M/s
201300200  88,839,464,508  		260.0 seconds elapsed 341.7M/s
505068400  241,847,443,968 		860.0 seconds elapsed 281.2M/s
1000000000 507,575,861,292
1006533300 511,171,793,022 		2180.0 seconds elapsed 234.5M/s

-- TBD
5000000000 2,891,512,530,078
5001011900 2,892,143,729,832		48265.0 seconds elapsed 59.9M/s
5174362000 3,000,375,710,610		50965.0 seconds elapsed 58.9M/s

*/


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
                    for (uint32_t i = f; i < sieve_length; i += f)
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
    for (; factor <= std::min<uint32_t>(sieve_length / 9, counts.size()); factor++) {
        if (factor <= MAX_BASE_F && sieve_length % factor == 0)
            continue;

        auto count = counts[factor];
        assert(count != -1ul);
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

    for (; factor <= std::min<uint32_t>(sieve_length / 2, counts.size()); factor++) {
        auto count = counts[factor];
        assert(count != -1ul);
        uint32_t index = count * factor - start;
        auto add = count + factor;

        for (; index < sieve_length; index += factor) {
            sums[index] += add++;
        }

        counts[factor] = add - factor;
    }

    // Handles factors that appear at least once (but possible more times)
    for (; factor <= std::min<uint32_t>(counts.size(), sieve_length); factor++) {
        auto count = counts[factor];
        assert(count != -1ul);
        uint32_t index = count * factor - start;
        auto add = count + factor;

        for (; index < sieve_length; index += factor) {
            sums[index] += add++;
        }

        counts[factor] = add - factor;
    }

    // Handles larger factors that can only appear once
    for (; factor <= counts.size(); factor++) {
        auto count = counts[factor];
        assert(count != -1ul);
        uint32_t index = count * factor - start;

        if (index < sieve_length) {
            // count = number / factor
            sums[index] += count + factor;
            counts[factor]++;
        }
    }

    start += sieve_length;
    return sums;
}

class A217259 {
    public:
        A217259(uint64_t start, uint64_t stop, uint64_t n)
            : START(start), STOP(stop), SEGMENT(n) {}

        bool test_match(uint16_t dist, uint64_t i);
        void print_match(uint16_t dist, uint64_t i);
        void print_match_and_test(uint16_t dist, uint64_t i) {
            print_match(dist, i);
            assert(test_match(dist, i));
        }

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

        const uint64_t START;
        const uint64_t STOP;
        const uint64_t SEGMENT;
        // DIST = 2 (A217259), 6 (A054903), 7 (A063680), 8 (A059118), 12 (A054902)
        const uint32_t MAX_DIST = 11;

        int64_t found = 0;
        int64_t found_composite = 0;
        int64_t print_mult = 1;
        uint64_t next_time = 5;

        std::chrono::time_point<std::chrono::system_clock> S =
            std::chrono::system_clock::now();

        // Guard for results
        std::mutex g_control;
        std::condition_variable work_ready;

        // Multithreaded results output
        vector<std::pair<uint64_t,std::unique_ptr<vector<uint64_t>>>> results;
};


bool A217259::test_match(uint16_t dist, uint64_t i) {
    // Already known terms
    // DIST = 2 (A217259), 6 (A054903), 7 (A063680), 8 (A059118), 12 (A054902)
    if (dist == 2)
        for (uint64_t known : {435, 8576, 8826})
            if ( (i+1) == known)
                return true;

    if (dist == 4)
        for (uint64_t known : {
                305635357,})
            if (i == known)
                return true;

    if (dist == 6)
        for (uint64_t known : std::initializer_list<uint64_t>{
                104, 147, 596, 1415, 4850, 5337, 370047, 1630622,
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
        if (mpz_probab_prime_p(t.get_mpz_t(), 20) != 2) {
            printf("\n\nNEW! %u, %lu | %lu not prime!\n\n", dist, i, i);
            return true;
        }
    }
    return true;
}

void A217259::print_match(uint16_t dist, uint64_t i) {
    found += 1;
    // So we can verify against A146214
    if ((found-found_composite) == 10*print_mult) {
        printf("\t%10ld'th prime pair: %'lu\n", found-found_composite, i);
        print_mult *= 10;
    } else if (found <= 10 || found % print_mult == 0) {
        printf("%-10ld %'-16lu\n", found, i);
    } else if (found % 100 == 0) {
        // Avoid calling sys_clock on every find.
        std::chrono::duration<double> elapsed = std::chrono::system_clock::now() - S;
        if (elapsed.count() > next_time) {
            // TODO add last interval rate
            float rate = (i - START) / elapsed.count() / 1e6;
            printf("%-10ld %'-16lu\t\t%.1f seconds elapsed %.1fM/s\n",
                    found, i, elapsed.count(), rate);
            // 5,10,15,95,100,110,120,130,...300,400,420,440
            next_time += 5 * (1 + (next_time >= 100)) * (1 + (next_time >= 300));
        }
    }
}


void A217259::iterate() {
    // sigma(start-), sigma(start-1)
    uint64_t last_sigmas[MAX_DIST];
    for (size_t i = 0; i < MAX_DIST; i++) {
        last_sigmas[i] = -1ul;
    }

    auto sieve = SegmentedSieveSigma(START, SEGMENT);

    for (uint64_t start = START; start <= STOP; start += SEGMENT) {
        //auto sigmas = SegmentedSieveOfSigma(start, SEGMENT);
        auto sigmas = sieve.next(start);

        // for (auto i = 0ul; i < 10; i++) {
        //     printf("\t%lu = %lu | %lu\n", i, i + start, sigmas[i]);

        for (size_t i = 0; i < MAX_DIST; i++) {
            for (size_t d = MAX_DIST - i; d <= MAX_DIST; i++) {
                // 9 is 1 distance from 0
                // 8 is 2 distance from 0
                assert(i + d >= MAX_DIST);

                if (last_sigmas[i] == sigmas[i + d - MAX_DIST]) {
                    if (last_sigmas[i] == 0)
                        print_match(d, start - 1);
                    else
                        print_match_and_test(d, start - 1);
                }
            }
        }

        for (uint32_t i = 1; i < SEGMENT-1; i++) {
            for (size_t d = 2; d <= MAX_DIST; i++) {
                if (sigmas[i] == sigmas[i + d]) {
                    if (sigmas[i] == 0 && sigmas[i + d] == 0)
                        print_match(d, start + i);
                    else {
                        printf("\tsigmas[%lu] = %lu, sigmas[%lu] = %lu\n",
                            start+i, sigmas[i], start+i+d, sigmas[i+d]);
                        print_match_and_test(d, start + i);
                    }
                }
            }
        }

        for (size_t i = 0; i < MAX_DIST; i++)
            last_sigmas[i] = sigmas[SEGMENT-MAX_DIST+i];
    }
}


void A217259::worker_thread() {
    #pragma omp parallel for schedule(dynamic, 4)
    for (uint64_t start = START; start <= STOP; start += (SEGMENT - MAX_DIST)) {
        // Calculate results
        vector<uint64_t> sigmas = SegmentedSieveOfSigma(start, SEGMENT);

        vector<uint64_t> terms;
        for (uint32_t i = 0; i < SEGMENT-MAX_DIST; i++) {
            for (uint16_t dist = 2; dist <= MAX_DIST; dist++) {
                if (sigmas[i] == sigmas[i + dist]) {
                    /**
                     * Old code verified i+1 and i-1 are prime explicitly.
                     * It's instant to instead check sigmas[i] == sigmas[j] == 0
                     * Can spot check with twin prime count (A146214)
                     */
                    if (sigmas[i] != 0 || sigmas[i+dist] != 0) {
                        printf("\tsigmas[%lu] = %lu, sigmas[%lu] = %lu\n",
                            start+i, sigmas[i], start+i+dist, sigmas[i+dist]);

                        assert(test_match(dist, start + i));
                    }

                    // TODO include dist
                    terms.push_back(start + i);
                }
            }
        }

        if ( (start+SEGMENT)  % 2'000'000'000l < SEGMENT ) {
            printf("\t\tComplete up to %'lu\n", start + SEGMENT - 1);
        }

        // Wait till I can safely queue results.
        std::unique_lock<std::mutex> guard(g_control);

        //printf("\t\tWork %lu -> %lu ready\n", start, terms.size());
        results.emplace_back(start, std::make_unique<vector<uint64_t>>(std::move(terms)));

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
    while (start <= STOP) {
        if (!guard) {
            guard.lock();
        }

        if (results.empty()) {
            // Release lock and wait for worker_thread
            work_ready.wait(guard);
            continue;
        }

        vector<uint64_t> *term_ptr = nullptr;

        // Check if start in results
        for (size_t i = 0; i < results.size(); i++) {
            if (results[i].first == start) {
                //printf("\tresults[%lu] = %lu\n", i, results[i].first);
                term_ptr = results[i].second.release();
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

        for (auto t : *term_ptr) {
            // TODO include dist
            print_match(2, t);
        }

        // Overlap by segments by MAX_DIST, because easier than saving final X values.
        start += SEGMENT - MAX_DIST;

        // Relock so that loop starts with lock
        guard.lock();
    }
}


int main() {
    // Allow comma seperators
    setlocale(LC_NUMERIC, "");

    printf("Compiled with GMP %d.%d.%d\n",
        __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL);

    uint64_t SEGMENT = 360360; //1 << 17;
    uint64_t START = 0; //(1'000'000'000 / SEGMENT + 1) * SEGMENT;
    uint64_t STOP = 1e13;

    A217259 runner(START, STOP, SEGMENT);

    // For single-threaded
    //runner.iterate();

    // For multi-threaded
    runner.multithreaded_iterate();
}
