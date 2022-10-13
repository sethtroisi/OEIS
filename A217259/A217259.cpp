// g++ -g -O3 --std=c++17 -Werror -Wall A217259.cpp -lgmpxx -lgmp -pthread -fopenmp && time ./a.out

/*

multi-threaded(5) 202
10000000   3,285,915,300
50000000   19,358,092,098
100000000  41,375,647,278
200000000  88,214,274,738
400000000  187,676,965,350
500000000  239,211,160,050
1000000000 507,575,861,292
2000000000 1,075,045,451,538
2000637100 1,075,414,804,482		10045.0 seconds elapsed 107.1M/s
           2,000,449,744,638		24025.0 seconds elapsed 83.3M/s
           3,749,496,713,070		59625.0 seconds elapsed 62.9M/s

10000000   3,285,915,300
100000000  41,375,647,278
500000000  239,211,160,050
1000000000 507,575,861,292
2000000000 1,075,045,451,538
2000419800 1,075,288,723,350		9505.0 seconds elapsed 113.1M/s

multi-threaded(4) 2022/10/13 results

10000000   3,285,915,300
13256600   4,488,240,300   		10.0 seconds elapsed 448.8M/s
100000000  41,375,647,278
103978300  43,179,449,502  		145.0 seconds elapsed 297.8M/s
200732200  88,566,518,868  		365.0 seconds elapsed 242.6M/s
304745000  139,593,089,568 		665.0 seconds elapsed 209.9M/s
400398500  187,880,326,968 		1005.0 seconds elapsed 186.9M/s
501047800  239,755,247,748 		1405.0 seconds elapsed 170.6M/s
	1000000000'th twin prime: 507,575,862,527
1002916700 509,181,043,518 		4005.0 seconds elapsed 127.1M/s
2002129700 1,076,282,265,762		11545.0 seconds elapsed 93.2M/s
3000000000 1,666,211,545,308
3000187400 1,666,323,730,872		21605.0 seconds elapsed 77.1M/s
3553850500 2,000,656,050,948		28085.0 seconds elapsed 71.2M/s
4000000000 2,273,005,950,738


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
 * Single threaded it would make sense to keep next_index array over isqrt?
 * O(read 4MB) < O(1 division + 1 multiplication)
 *
 * For very small factors (2...50)
 *   If SEGMENT was a multiple of lcm(2,3,4,5,6,7), could have a base + increment value
 *   [60,66) init
 *        [2+60/2 + 3+60/3, 0, 2+62/2, 3+63/3, 2+64/2, 0) =
 *        [32+23, 0, 33, 24, 34, 0)
 *        [55, 0, 33, 24, 34, 0)
 *   [66,72) init
 *        [2+66/2 + 3+66/3, 0, 2+68/2, 3+69/3, 2+70/2, 0) =
 *        [35+25, 0, 36, 26, 37, 0)
 *        [60, 0, 36, 26, 37, 0)
 *            delta from first [5, 0, 3, 2, 3, 0]
 *                            =[3, 0, 3, 0, 3, 0] = 6/2 every 2nd
 *                            +[2, 0, 0, 2, 0, 0] = 6/3 every 3rd
 *   [72,78) init
 *        [65, 0, 39, 28, 40, 0)
 *
 *   lcm(2, ..., 12) = 27720  = 14.8 bits
 *   lcm(2, ..., 15) = 360360 = 18.5 bits
 *   1/2 + 1/3 + 1/4 + 1/5 ... 1/15 = 2.318
 *
 *   Can reduce 2.318 * N random writes to 2*N (update, copy) linear writes
 *        It's kinda more like 2.318 * 2 because of index?
 *   On the order of 7 * N total writes so this could be a big win?
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
            // count = number / factor
            sums[index           ]  += add++;
            sums[index +   factor]  += add++;
            sums[index + 2*factor]  += add++;
            sums[index + 3*factor]  += add++;
            index += factor<<2;
        }

        for (; index < N; index += factor) {
            sums[index] += add++;
        }
    }

    // Handles factors that can appear more than once
    for (; factor <= std::min(isqrt, N); factor++) {
        uint64_t count = start_m_1 / factor + 1;
        uint32_t index = count * factor - start;
        uint64_t add = factor + count;

        for (; index < N; index += factor) {
            sums[index] += add++;
        }
    }

    // Handles larger factors which can only appear once
    for (; factor <= isqrt; factor++) {
        uint64_t count = start_m_1 / factor + 1;
        uint32_t index = count * factor - start;
        uint64_t add = factor + count;

        if (index < N) {
          // count = number / factor
          sums[index] += add;
        }
    }

    return sums;
}

class SegmentedSieveSigma {
    public:
        SegmentedSieveSigma(uint64_t start, uint64_t length) : sieve_length(length) {
            sums.resize(sieve_length);
            if (start > 0) {
                jump_to(start);
            }
        }

        void jump_to(uint64_t new_start);
        const vector<uint64_t>& next(uint64_t verify_start);

    private:
        const uint32_t sieve_length;

        // Start of the next interval
        uint64_t start = 0;

        // last number in offsets array, isqrt = floor(sqrt(start-1));
        uint32_t isqrt = 1;
        vector<std::pair<uint32_t,uint64_t>> offsets = {{0,0}, {0,0}};

        vector<uint64_t> sums;
};

void SegmentedSieveSigma::jump_to(uint64_t new_start) {
    assert(new_start > start);

    // Add factors, f, while f*f < start
    auto factors = calc_isqrt(new_start - 1);
    assert((factors+1) >= offsets.size());

    // Possible add new empty offsets
    offsets.resize(factors + 1);

    for (uint64_t f = 2; f <= factors; f++) {
        // -new_start % f
        uint64_t count = (new_start-1) / f + 1;
        uint32_t index = count * f - new_start;
        offsets[f] = {index, f + count};
    }
 }

const vector<uint64_t>& SegmentedSieveSigma::next(uint64_t verify_start) {
    assert(start == verify_start);

    //auto sums = vector<uint64_t>(sieve_length, 0);
    std::memset(&sums[0], 0, sieve_length * sizeof(sums[0]));

    auto past = start + sieve_length;

    // For all new factors
    for(uint64_t f = isqrt+1, f2 = f * f; f2 < past; f++, f2 = f*f) {
        assert(f2 >= start);

        uint32_t index = f2 - start;
        assert(index < sieve_length);

        // n = factor^2 only adds factor because n / factor = factor.
        sums[index] += f;

        // Add offset for lower loop to handle
        assert(offsets.size() == f);

        isqrt = f;
        offsets.push_back({index + f, f + (f + 1)});
    }
    assert((uint64_t) (isqrt+1)*(isqrt+1) >= past);

    // Hand unrolled loops for very small factor
    uint64_t factor = 2;
    for (; factor <= std::min(isqrt, sieve_length/5); factor++) {
        auto [index, add] = offsets[factor];

        // Loop unrolled 4x for small factors
        for (; index + (factor<<2) < sieve_length; ) {
            // count = number / factor
            sums[index           ]  += add++;
            sums[index +   factor]  += add++;
            sums[index + 2*factor]  += add++;
            sums[index + 3*factor]  += add++;
            index += factor<<2;
        }

        for (; index < sieve_length; index += factor) {
            sums[index] += add++;
        }

        offsets[factor] = {index - sieve_length, add};
    }

    // Handles factors that appear at least once (but possible more times)
    for (; factor <= std::min(isqrt, sieve_length); factor++) {
        auto [index, add] = offsets[factor];

        for (; index < sieve_length; index += factor) {
            sums[index] += add++;
        }

        offsets[factor] = {index - sieve_length, add};
    }

    // Handles larger factors that can only appear once
    for (; factor <= isqrt; factor++) {
        auto [index, add] = offsets[factor];

        if (index < sieve_length) {
          // count = number / factor
          sums[index] += add;
          offsets[factor] = {index + factor, add + 1};
        } else {
          offsets[factor] = {index - sieve_length, add};
        }
    }

    start += sieve_length;
    return sums;
}

class A217259 {
    public:
        A217259(uint64_t start, uint64_t stop, uint64_t n)
            : START(start), STOP(stop), SEGMENT(n) {}

        bool test_match(uint64_t mid);
        void print_match(uint64_t mid);
        void print_match_and_test(uint64_t mid) {
            print_match(mid);
            assert(test_match(mid));
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

        int64_t found = 0;
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


bool A217259::test_match(uint64_t mid) {
    // Special already know terms
    if (mid == 435 || mid == 8576 || mid == 8826)
        return true;

    mpz_class mid_m_1 = mid - 1;
    if (mpz_probab_prime_p(mid_m_1.get_mpz_t(), 20) != 2) {
        printf("NEW! %lu | %lu - 1 not prime!\n", mid, mid + 1);
        return false;
    }
    mpz_class mid_p_1 = mid + 1;
    if (mpz_probab_prime_p(mid_p_1.get_mpz_t(), 20) != 2) {
        printf("NEW! %lu | %lu + 1 not prime!\n", mid, mid + 1);
        return false;
    }
    return true;
}

void A217259::print_match(uint64_t mid) {
    found += 1;
    // So we can verify against A146214
    if ((found-3) == 10*print_mult) {
        if (mid > 8826)
          printf("\t%10ld'th twin prime: %'lu\n", found-3, mid-1);
        print_mult *= 10;
    } else if (found <= 10 || found % print_mult == 0) {
        printf("%-10ld %'-16lu\n", found, mid);
    } else if (found % 100 == 0) {
        // Avoid calling sys_clock on every find.
        std::chrono::duration<double> elapsed = std::chrono::system_clock::now() - S;
        if (elapsed.count() > next_time) {
            // TODO add last interval rate
            float rate = (mid - START) / elapsed.count() / 1e6;
            printf("%-10ld %'-16lu\t\t%.1f seconds elapsed %.1fM/s\n",
                    found, mid, elapsed.count(), rate);
            // 5,10,15,95,100,110,120,130,...300,400,420,440
            next_time += 5 * (1 + (next_time >= 100)) * (1 + (next_time >= 300));
        }
    }
}


void A217259::iterate() {
    // sigma(start-2), sigma(start-1)
    uint64_t last_sigmas[2] = {-1ul, -1ul};

    if (START <= 1)
        printf("\tCAN'T CHECK %lu or %lu\n", START, START+1);

    auto sieve = SegmentedSieveSigma(START, SEGMENT);

    for (uint64_t start = START; start <= STOP; start += SEGMENT) {
        //auto sigmas = SegmentedSieveOfSigma(start, SEGMENT);
        auto sigmas = sieve.next(start);

        if (sigmas[0] == last_sigmas[0]) {
            if (sigmas[0] == 0)
                print_match(start - 1);
            else
                print_match_and_test(start - 1);
        }

        if (sigmas[1] == last_sigmas[1]) {
            if (sigmas[1] == 0)
                print_match(start);
            else
                print_match_and_test(start);
        }

        for (uint32_t i = 1; i < SEGMENT-1; i++) {
            if (sigmas[i+1] == sigmas[i-1]) {
                if (sigmas[i+1] == 0 && sigmas[i-1] == 0)
                    print_match(start + i);
                else
                    print_match_and_test(start + i);
            }
        }

        last_sigmas[0] = sigmas[SEGMENT-2];
        last_sigmas[1] = sigmas[SEGMENT-1];
    }
}


void A217259::worker_thread() {
    #pragma omp parallel for schedule(dynamic, 4)
    for (uint64_t start = START; start <= STOP; start += (SEGMENT - 2)) {
        // Calculate results
        vector<uint64_t> sigmas = SegmentedSieveOfSigma(start, SEGMENT);

        vector<uint64_t> terms;
        for (uint32_t i = 1; i < SEGMENT-1; i++) {
            if (sigmas[i+1] == sigmas[i-1]) {
                /**
                 * Used to verify i+1 and i-1 are prime explicitly.
                 * Faster to check sigmas[i+1] == sigmas[i-1] == 0
                 * Then spot check with twin prime count
                 */
                if (sigmas[i-1] != 0 || sigmas[i+1] != 0) {
                  assert(test_match(start + i));
                }
                terms.push_back(start + i);
            }
        }

        //if ( (start+SEGMENT)  % 2'000'000'000l < SEGMENT ) {
        //    printf("\t\t[%lu, %lu) complete\n", start, start + SEGMENT);
        //}

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
    // Overlap by segments by two, because easier
    // keep queue of open intervals,

    if (START <= 1)
        printf("\tCAN'T CHECK %lu or %lu\n", START, START+1);

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
            print_match(t);
        }

        start += SEGMENT - 2;

        // Relock so that loop starts with lock
        guard.lock();
    }
}


int main() {
    // Allow comma seperators
    setlocale(LC_NUMERIC, "");

    printf("Compiled with GMP %d.%d.%d\n",
        __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL);

    uint64_t START = 0;
    uint64_t SEGMENT = 1 << 17;
    uint64_t STOP = 1e13;

    A217259 runner(START, STOP, SEGMENT);

    // For single-threaded
    runner.iterate();

    // For multi-threaded
    //runner.multithreaded_iterate();
}
