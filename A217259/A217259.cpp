// g++ -g --std=c++14 -O3 -Werror -Wall A217259.cpp -lgmpxx -lgmp -pthread -fopenmp && time ./a.out

/*

single-threaded
100000000 41,375,647,278
100050453 41,398,462,908		1350.0 seconds elapsed 30.7M/s
300000000 137,227,128,108       ~5130
400000000 187,676,965,350       ~7415
424182497 200,050,245,012       7995.0 seconds elapsed 25.0M/s

multi-threaded (5 threads)
10000000 3,285,915,300
10174200 3,349,511,130		    20.0 seconds elapsed 167.5M/s
50000000 19,358,092,098
50395900 19,526,164,830		    125.0 seconds elapsed 156.2M/s
100000000 41,375,647,278
100028700 41,388,642,768		280.0 seconds elapsed 147.8M/s

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
 * Calculate sigma[i] for i in range(start, start+n).
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

    // Adjust to include n as a divisor of n
    for (uint64_t i = 0; i < N; i++) {
        sigmas[i] += start + i;
    }

    // isqrt with gmp
    mpz_class sqrt = past-1;
    mpz_sqrt(sqrt.get_mpz_t(), sqrt.get_mpz_t());
    uint64_t isqrt = mpz_get_ui(sqrt.get_mpz_t());

    assert( isqrt * isqrt < past );
    assert( (isqrt+1) * (isqrt+1) >= past );

    for (uint64_t factor = 2; factor <= isqrt; factor++) {
        auto f2 = factor * factor;
        assert(f2 < past);

        uint64_t ceil, next_index;

        if (f2 >= start) {
            assert(f2 - start < N);
            // n=factor^2 only gets to count factor, not factor + n/factor
            sigmas[f2 - start] += factor;

            ceil = factor + 1;
            next_index = f2 + factor - start;
        } else {
            // ceil(start / factor)
            ceil = (start-1) / factor + 1;
            next_index = ceil * factor - start;
        }

        uint64_t count = ceil;
        for (uint64_t index = next_index; index < N; index += factor, count++) {
            // count = number / factor
            sigmas[index] += factor + count;
        }
    }

    return sigmas;
}


bool test_match(uint64_t mid) {
    if (mid == 435 or mid == 8576 or mid == 8826)
        return true;

    // Verify mid-1 and mid+1 are prime
    mpz_class mid_m_1 = mid - 1;
    mpz_class mid_p_1 = mid + 1;
    if((mpz_probab_prime_p(mid_m_1.get_mpz_t(), 20) != 2) ||
       (mpz_probab_prime_p(mid_p_1.get_mpz_t(), 20) != 2)) {
        printf("NEW! %lu | %lu or %lu\n", mid, mid - 1, mid + 1);
        exit(1);
    }
    return true;
}

void print_match(uint64_t mid) {
    static uint32_t found = 0;
    static uint64_t print_mult = 1;
    static auto S = std::chrono::system_clock::now();
    static auto next_time = 5;


    // TODO commas in mid

    found += 1;
    if (found % print_mult == 0) {
        printf("%-10d %'-16lu\n", found, mid);
        if (found == 10 * print_mult)
            print_mult *= 10;
    } else if (found % 100 == 0) {
        // Avoid calling sys_clock on every find.
        std::chrono::duration<double> elapsed = std::chrono::system_clock::now() - S;
        if (elapsed.count() > next_time) {
            float rate = mid / elapsed.count() / 1e6;
            printf("%-10d %'-16lu\t\t%.1f seconds elapsed %.1fM/s\n",
                    found, mid, elapsed.count(), rate);
            next_time += 5;
        }
    }
}

void print_match_and_test(uint64_t mid) {
    print_match(mid);
    test_match(mid);
}

void iterate(uint64_t START, uint64_t STOP, uint64_t SEGMENT) {
    // sigma(start-2), sigma(start-1)
    uint64_t last_sigmas[2] = {0, 0};

    if (START > 0)
        printf("\tCAN'T CHECK %lu or %lu\n", START, START+1);

    for (uint64_t start = START; start <= STOP; start += SEGMENT) {
        auto sigmas = SegmentedSieveOfSigma(start, SEGMENT);

        if (sigmas[0] - last_sigmas[0] == 2)
            print_match_and_test(start - 1);

        if (sigmas[1] - last_sigmas[1] == 2)
            print_match_and_test(start);

        for (uint32_t i = 1; i < SEGMENT-1; i++) {
            if (sigmas[i+1] - sigmas[i-1] == 2) {
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


void worker_thread(uint64_t START, uint64_t STOP, uint64_t SEGMENT) {
    // pass handles to vector around?

    #pragma omp parallel for schedule(dynamic, 1)
    for (uint64_t start = START; start <= STOP; start += (SEGMENT - 2)) {
        // Calculate results
        vector<uint64_t> sigmas = SegmentedSieveOfSigma(start, SEGMENT);

        vector<uint64_t> terms;
        for (uint32_t i = 1; i < SEGMENT-1; i++) {
            if (sigmas[i+1] - sigmas[i-1] == 2) {
                assert(test_match(start + i));
                terms.push_back(start + i);
            }
        }

        // Wait till I can queue results
        std::unique_lock<std::mutex> guard(g_control);

        //printf("\t\tWork %lu -> %lu ready\n", start, terms.size());
        results.emplace_back(start, std::make_unique<vector<uint64_t>>(std::move(terms)));

        // Let iteratator advance
        guard.unlock();
        work_ready.notify_one();
        //printf("\t\tresults %.4f, guard held for %.4f\n",
        //        calculation.count(), guarded.count());
    }
    printf("\t\tAll work finished!\n");
}

void multithreaded_iterate(uint64_t START, uint64_t STOP, uint64_t SEGMENT) {
    // Overlap by segments by two, because easier
    // keep queue of open intervals,

    if (START > 0)
        printf("\tCAN'T CHECK %lu or %lu\n", START, START+1);

    // Wait for added item
    std::unique_lock<std::mutex> guard(g_control);
    // Wait for some work ready
    work_ready.wait(guard);
    printf("\tSome work ready!\n");

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

        //printf("\tlooking for work(%lu)\n", results.size());
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
        //if (!results.empty())
        //    printf("\tNow %lu results\n", results.size());

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

        // Relock so that look starts with lock
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
    uint64_t STOP = 200e9;
    //iterate(START, STOP, SEGMENT);

    std::thread t1(worker_thread, START, STOP, SEGMENT);
    std::thread t2(multithreaded_iterate, START, STOP, SEGMENT);

    t1.join();
    t2.join();
}
