#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <numeric>
#include <ranges>
#include <string>
#include <utility>
#include <vector>

#include <gmpxx.h>
#include <omp.h>

//#include "cgbn_compute_cf.h"
#include "int128_compute_cf.h"

#ifndef _OPENMP
    // Fakes in case -fopenmp isn't passed
    int omp_get_max_threads() { return 1; }
    int omp_get_thread_num() { return 0; }
#endif

// Relates to fibonacci and max solution
size_t MAX_CF = 1'000;

using std::chrono::duration;
using std::chrono::high_resolution_clock;
using std::pair;
using std::vector;

vector<uint32_t> get_primes(uint32_t n) {
    vector<uint32_t> primes;

    for (uint32_t i = 2; i <= n; i++) {
        bool is_prime = true;
        for (uint32_t j = 2; j * j <= i; ++j) {
            if (i % j == 0) {
                is_prime = false;
                break;
            }
        }

        if (is_prime) {
            primes.push_back(i);
        }
    }
    return primes;
}


vector<mpz_class> power_set(const vector<uint32_t>& set) {
    size_t n = set.size();
    assert( n < 30 );
    size_t size = 1 << n; // Equivalent to 2^n
    // Set all entries to 1.
    vector<mpz_class> ps;
    ps.reserve(size);
    ps.push_back(1); // For entry 0

    for (size_t i = 1; i < size; i++) {
        // Guaranteed to have some bits set
        size_t j = __builtin_ctzll(i);
        assert( i & (1 << j) );
        mpz_class p = set[j] * ps[i ^ (1 << j)];
        ps.push_back(p);
    }

    std::sort(ps.begin(), ps.end());
    // Do the hardest part first!
    std::reverse(ps.begin(), ps.end());
    return ps;
}


__uint128_t from_mpz_class(const mpz_class& t) {
    // Awkward hack two limb x into t
    return (((__uint128_t) mpz_getlimbn(t.get_mpz_t(), 1)) << 64) | mpz_getlimbn(t.get_mpz_t(), 0);
}


// https://mathworld.wolfram.com/PeriodicContinuedFraction.html#eqn2
// for squarefree D, 0 < ak < 2 * sqrt(n)
pair<bool, uint64_t> continued_fraction_sqrt_126_pessemistic(mpz_class x_in) {
    assert( mpz_sizeinbase(x_in.get_mpz_t(), 2) <= 126 );

    mpz_class t = sqrt(x_in);
    __uint128_t a0 = mpz_get_ui(t.get_mpz_t());
    __uint128_t x = from_mpz_class(x_in);
    __uint128_t b = a0;
    __uint128_t c = x - b*b;
    __uint128_t a = (a0 << 1) / c;

    uint64_t i = 2; // a0, a

    __uint128_t two_a0 = a0 << 1;
    for (; i <= MAX_CF && a != two_a0; ) {
        b = a*c - b;
        c = (x - b*b) / c;
        a = (a0 + b) / c;

        // 1 <= b <= a0
        // 1 <= c <= a0 + b
        // c | (x - b*b)
        ++i;
    }
    return {a == two_a0, i};
}


/**
 * Given prime p, P with p <= P
 * handle Q = powerset([2, 3, 7, 11, p])
 */
void StormersTheorem(uint32_t p, uint32_t P) {
    auto primes = get_primes(P);
    assert( P == primes.back() );

    vector<int> p_index(P+1, 0);
    for (size_t i = 0; i < primes.size(); i++) p_index[primes[i]] = i;

    int p_i = p_index[p];
    assert( primes[p_i] == p );

    // Q, CF < MAX_CF, GPU version
    static uint64_t total[3] = {};
    double times[2] = {};

    /**
     * Minimize memory usage by breaking Q' into half
     * LOWER size gives us update frequency for fancy printing
     * UPPER size helps with multithreading.
     */
    const int32_t HIGH_PRIMES = p_i < 15 ? p_i : 4 * p_i / 5;
    assert( 0 <= HIGH_PRIMES && HIGH_PRIMES <= p_i);
    vector<uint32_t> primes_low(primes.begin(), primes.begin() + (p_i - HIGH_PRIMES));
    vector<uint32_t> primes_high(primes.begin() + (p_i - HIGH_PRIMES), primes.begin() + p_i);
    const vector<mpz_class> Q_low = power_set(primes_low);
    const vector<mpz_class> Q_high = power_set(primes_high);

    vector<pair<__uint128_t, __uint128_t>> temp_Q(Q_high.size(), {0, 0});
    vector<uint32_t> valid(Q_high.size(), 0);

    // Inner loop temporaries
    mpz_class D, q, x_1, y_1, x_n, y_n, x_np1, y_np1, x, y, t;
    //CgbnPessemisticCf gpu_tester(Q_high.size());
    PessemisticCf gpu_tester(Q_high.size());


    printf("p: %u\n", p);
    for (mpz_class Q_1 : Q_low) {
        // Always include p as a convince multiply into Q_1 here
        // This means q=1 is skipped but that's fine as it doesn't generate solutions.
        Q_1 *= p;

        vector<uint64_t> local_counts(omp_get_max_threads(), 0);

        auto gpu_start = high_resolution_clock::now();
        #pragma omp parallel for schedule(dynamic) private(D, t)
        for (size_t i = 0; i < Q_high.size(); i++) {
            const mpz_class& Q_2 = Q_high[i];
            D = Q_1 * Q_2;
            t = sqrt(D);
            temp_Q[i] = {from_mpz_class(D), from_mpz_class(t)};
        }

        gpu_tester.run(MAX_CF, temp_Q, valid, false);
        for (size_t i = 0; i < Q_high.size(); i++) {
            if (valid[i] > 0)
                total[2]++;
        }
        duration<double> gpu_time = high_resolution_clock::now() - gpu_start;

        auto cpu_start = high_resolution_clock::now();
        /*
        #pragma omp parallel for schedule(dynamic) \
            private(D, q, x_1, y_1, x_n, y_n, x_np1, y_np1, x, y)
        //for (const mpz_class& Q_2 : Q_high) {
        for (size_t i = 0; i < Q_high.size(); i++) {
            const mpz_class& Q_2 = Q_high[i];
            q = Q_1 * Q_2;

            // Lucas computes D = q which generates other interesting numbers
            // Lehmer used D = 2 * q which only generates A002071
            D = q; // * 2;

            uint64_t &count = local_counts[omp_get_thread_num()];
            bool is_small = (mpz_sizeinbase(D.get_mpz_t(), 2) <= 126);
            assert( is_small );

            auto t = continued_fraction_sqrt_126_pessemistic(D);
            if (t.first) {
                count++;
            }

            //gmp_printf("%Zd -> {%d, %lu} vs {%lu}\n", D, t.first, t.second, valid[i]);
        }
        //*/
        duration<double> cpu_time = high_resolution_clock::now() - gpu_start;

        total[0] += Q_high.size();
        for (int i = 0; i < omp_get_max_threads(); i++) {
            total[1] += local_counts[i];
        }
        times[0] += cpu_time.count();
        times[1] += gpu_time.count();

        printf("\t%lu -> %lu vs %lu (%.1f vs %.1f)\n",
                total[0], total[1], total[2], times[0], times[1]);
    }
}


int main(int argc, char** argv) {
    assert(argc == 3);
    assert(mp_bits_per_limb == 64);

    int n = atol(argv[1]);
    MAX_CF = atol(argv[2]);

    auto primes = get_primes(n);
    auto P = primes.back();
    if (n < 0 || ((unsigned) n != P)) {
        printf("Usage: %s P CF\n", argv[0]);
        printf("P(%d) must be prime\n", n);
        exit(1);
    }

    if (MAX_CF < 10 || MAX_CF >= 1'000'000'000) {
        printf("Usage: %s P CF\n", argv[0]);
        printf("CF(%lu) must be >= 10 and < 1B\n", MAX_CF);
        exit(1);
    }


    if (1) {
        double primorial_P = std::accumulate(primes.begin(), primes.end(), 1.0, std::multiplies<>{});
        auto d = log2(primorial_P);

        if (1) { // Log of product of primes, tells us about square(2*q)
            printf("|Primes %d ... %d| = %lu -> log2(P) = %.1f | MAX_CF = %lu\n",
                   primes.front(), primes.back(), primes.size(), d, MAX_CF);
        }
    }

    for (uint32_t p_i = 0; p_i < primes.size(); p_i++) {
        auto p = primes[p_i];
        StormersTheorem(p, P);
    }
}
