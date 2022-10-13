// g++ --std=c++14 -O3 -Werror -Wall A217259.cpp -lgmpxx -lgmp -pthread -fopenmp && time ./a.out

/*

single-threaded (testing both mid-1 and mid+1)
100000000 41,375,647,278
100050453 41,398,462,908		1350.0 seconds elapsed 30.7M/s
300000000 137,227,128,108       ~5130
400000000 187,676,965,350       ~7415
424182497 200,050,245,012       7995.0 seconds elapsed 25.0M/s

multi-threaded (5 threads) (testing both)
10000000 3,285,915,300
10174200 3,349,511,130		    20.0 seconds elapsed 167.5M/s
50000000 19,358,092,098
50395900 19,526,164,830		    125.0 seconds elapsed 156.2M/s
100000000 41,375,647,278
100028700 41,388,642,768		280.0 seconds elapsed 147.8M/s

multi-threaded (5 threads)
10000000   3,285,915,300
16987200   5,901,309,372   		15.0 seconds elapsed 393.4M/s
50000000   19,358,092,098
54779200   21,396,826,278  		65.0 seconds elapsed 329.2M/s
100000000  41,375,647,278
100475800  41,591,298,180  		145.0 seconds elapsed 286.8M/s
200000000  88,214,274,738
203188600  89,748,612,852  		375.0 seconds elapsed 239.3M/s
400000000  187,676,965,350
400595700  187,981,263,330 		965.0 seconds elapsed 194.8M/s
500000000  239,211,160,050
502046500  240,274,643,628 		1325.0 seconds elapsed 181.3M/s
1000000000 507,575,861,292
1000393400 507,792,465,780 		3585.0 seconds elapsed 141.6M/s
2000000000 1,075,045,451,538
2000637100 1,075,414,804,482		10045.0 seconds elapsed 107.1M/s
           2,000,449,744,638		24025.0 seconds elapsed 83.3M/s
           3,749,496,713,070		59625.0 seconds elapsed 62.9M/s


multi-threaded (5 threads) - newish, no primality code
10000000   3,285,915,300
18966500   6,664,539,282   		15.0 seconds elapsed 444.3M/s
100000000  41,375,647,278
105196300  43,732,190,598  		135.0 seconds elapsed 323.9M/s
203639500  89,965,823,148  		335.0 seconds elapsed 268.6M/s
302426200  138,436,564,548 		585.0 seconds elapsed 236.6M/s
403039100  189,229,712,718 		885.0 seconds elapsed 213.8M/s
500000000  239,211,160,050
500281100  239,357,154,612 		1215.0 seconds elapsed 197.0M/s
1000000000 507,575,861,292
1001622100 508,468,670,418 		3415.0 seconds elapsed 148.9M/s
2000000000 1,075,045,451,538
2000419800 1,075,288,723,350		9505.0 seconds elapsed 113.1M/s

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

/**
 * Calculate ~sigma[i]~ aliquot_sum[i] for i in range(start, start+n).
 */
vector<uint64_t> SegmentedSieveOfSigma(uint64_t start, uint64_t N) {
    auto past = start + N;
    auto sigmas = vector<uint64_t>(N, 1);

    // HACK for sigmas(0) and sigmas(1)
    for (int i = 0; i < 2; i++) {
        if (start + i < 2) {
            sigmas[i] = 0;
        }
    }

    // Compute aliquot sum (exclude N) -> requires removing +2 below
    /*
    // Adjust to include n as a divisor of n
    for (uint64_t i = 0; i < N; i++) {
        sigmas[i] += start + i;
    }
    */

    // isqrt with gmp
    mpz_class sqrt = past-1;
    mpz_sqrt(sqrt.get_mpz_t(), sqrt.get_mpz_t());
    uint64_t isqrt = mpz_get_ui(sqrt.get_mpz_t());

    assert( isqrt * isqrt < past );
    assert( (isqrt+1) * (isqrt+1) >= past );

    /**
     * Handle the few factors with start < f^2 < start+N seperatly
     */
    for (; isqrt >= 2; isqrt--) {
        uint32_t factor = isqrt;
        uint64_t f2 = (uint64_t) factor * factor;
        if (f2 < start)
          break;

        assert(f2 < start + N);

        uint32_t index = f2 - start;
        // n = factor^2 only adds factor because n / factor = factor.
        sigmas[index] += factor;

        // n = factor * (factor+1)
        uint32_t aliquot_add = (factor) + (factor + 1);
        for (index += factor; index < N; index += factor) {
            sigmas[index] += aliquot_add++;
        }
    }

    /**
     * Single threaded it would make sense to keep next_index array over isqrt?
     * O(read 4MB) < O(1 division + 1 multiplication)
     */

    /**
     * NOTE!
     * Roughly half of all time is spent in very small factor loop
     *
     * If SEGMENT was a multiple of lcm(2,3,4,5,6,7), could have some mask add
     * (2+n/2 + 3+n/3 +4+n/4),0,(2+(n)/2+1),(3+n/3),(2+(n)/2+2 + 4+n/4+1
     * (two + three + four), 0, (++two), (++three), (++two + ++four)
     * [60,66) init
     *      [2+60/2 + 3+60/3, 0, 2+62/2, 3+63/3, 2+64/2, 0) =
     *      [32+23, 0, 33, 24, 34, 0)
     *      [55, 0, 33, 24, 34, 0)
     * [66,72) init
     *      [2+66/2 + 3+66/3, 0, 2+68/2, 3+69/3, 2+70/2, 0) =
     *      [35+25, 0, 36, 26, 37, 0)
     *      [60, 0, 36, 26, 37, 0)
     *          delta from first [5, 0, 3, 2, 3, 0]
     *                          =[3, 0, 3, 0, 3, 0] = 6/2 every 2nd
     *                          +[2, 0, 0, 2, 0, 0] = 6/3 every 3rd
     * [72,78) init
     *      [70, 0, 39, 28, 40, 0)
     *
     * lcm(2, ..., 12) = 27720  = 14.8 bits
     * lcm(2, ..., 15) = 360360 = 18.5 bits
     * 1/2 + 1/3 + 1/4 + 1/5 ... 1/15 = 2.318
     *
     * Can reduce 2.318 * N random writes to 2*N (update, copy) linear writes
     */
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
            sigmas[index           ]  += add++;
            sigmas[index +   factor]  += add++;
            sigmas[index + 2*factor]  += add++;
            sigmas[index + 3*factor]  += add++;
            index += factor<<2;
        }

        for (; index < N; index += factor) {
            sigmas[index] += add++;
        }
    }

    // Handles factors that can appear more than once
    for (; factor <= std::min(isqrt, N); factor++) {
        uint64_t count = start_m_1 / factor + 1;
        uint32_t index = count * factor - start;
        uint64_t add = factor + count;

        for (; index < N; index += factor) {
            sigmas[index] += add++;
        }
    }

    // Handles larger factors which can only appear once
    for (; factor <= isqrt; factor++) {
        uint64_t count = start_m_1 / factor + 1;
        uint32_t index = count * factor - start;
        uint64_t add = factor + count;

        if (index < N) {
          // count = number / factor
          sigmas[index] += add;
        }
    }

    return sigmas;
}


bool test_match(uint64_t mid) {
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

void print_match(uint64_t mid) {
    static int64_t found = 0;
    static int64_t print_mult = 1;
    static auto S = std::chrono::system_clock::now();
    static auto next_time = 5;
    static auto first_mid = mid;

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
            float rate = (mid - first_mid) / elapsed.count() / 1e6;
            printf("%-10ld %'-16lu\t\t%.1f seconds elapsed %.1fM/s\n",
                    found, mid, elapsed.count(), rate);
            // 5,10,15,95,100,110,120,130,...300,400,420,440
            next_time += 5 * (1 + (next_time > 100)) * (1 + (next_time > 300));
        }
    }
}

void print_match_and_test(uint64_t mid) {
    print_match(mid);
    assert(test_match(mid));
}

void iterate(uint64_t START, uint64_t STOP, uint64_t SEGMENT) {
    // sigma(start-2), sigma(start-1)
    uint64_t last_sigmas[2] = {-1ul, -1ul};

    if (START <= 1)
        printf("\tCAN'T CHECK %lu or %lu\n", START, START+1);

    for (uint64_t start = START; start <= STOP; start += SEGMENT) {
        auto sigmas = SegmentedSieveOfSigma(start, SEGMENT);

        if (sigmas[0] == last_sigmas[0]) {
            if (sigmas[0] == 1)
                print_match(start - 1);
            else
                print_match_and_test(start - 1);
        }

        if (sigmas[1] == last_sigmas[1]) {
            if (sigmas[1] == 1)
                print_match(start);
            else
                print_match_and_test(start);
        }

        for (uint32_t i = 1; i < SEGMENT-1; i++) {
            if (sigmas[i+1] == sigmas[i-1]) {
                if (sigmas[i+1] == 1 && sigmas[i-1] == 1)
                    print_match(start + i);
                else
                    print_match_and_test(start + i);
            }
        }

        last_sigmas[0] = sigmas[SEGMENT-2];
        last_sigmas[1] = sigmas[SEGMENT-1];
    }
}


// Guard for results
std::mutex g_control;
std::condition_variable work_ready;
vector<std::pair<uint64_t,std::unique_ptr<vector<uint64_t>>>> results;


// It's okay to ignore openmp pragma if not compiled with openmp.

void worker_thread(uint64_t START, uint64_t STOP, uint64_t SEGMENT) {
    #pragma omp parallel for schedule(dynamic, 4)
    for (uint64_t start = START; start <= STOP; start += (SEGMENT - 2)) {
        // Calculate results
        vector<uint64_t> sigmas = SegmentedSieveOfSigma(start, SEGMENT);

        vector<uint64_t> terms;
        for (uint32_t i = 1; i < SEGMENT-1; i++) {
            if (sigmas[i+1] == sigmas[i-1]) {
                /**
                 * Used to verify i+1 and i-1 are prime explicitly.
                 * Faster to check sigmas[i+1] == sigmas[i-1] == 1
                 * Then spot check with twin prime count
                 */
                if (sigmas[i-1] != 1 || sigmas[i+1] != 1) {
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

void multithreaded_iterate(uint64_t START, uint64_t STOP, uint64_t SEGMENT) {
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
    //iterate(START, STOP, SEGMENT);

    std::thread t1(worker_thread, START, STOP, SEGMENT);
    std::thread t2(multithreaded_iterate, START, STOP, SEGMENT);

    t1.join();
    t2.join();
}
