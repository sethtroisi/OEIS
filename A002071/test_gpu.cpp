#include <algorithm>
#include <cassert>
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

#include "cgbn_compute_cf.h"

#ifndef _OPENMP
    // Fakes in case -fopenmp isn't passed
    int omp_get_max_threads() { return 1; }
    int omp_get_thread_num() { return 0; }
#endif

// Relates to fibonacci and max solution
size_t MAX_CF = 1'000;

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
    uint64_t total[3];

    /**
     * Minimize memory usage by breaking Q' into half
     * LOWER size gives us update frequency for fancy printing
     * UPPER size helps with multithreading.
     */
    const int32_t LOW_PRIMES = p_i < 20 ? std::min<int32_t>(p_i, 2) : (p_i - 10) / 2;
    assert( 0 <= LOW_PRIMES && LOW_PRIMES <= p_i);
    vector<uint32_t> primes_low(primes.begin(), primes.begin() + LOW_PRIMES);
    vector<uint32_t> primes_high(primes.begin() + LOW_PRIMES, primes.begin() + p_i);
    const vector<mpz_class> Q_low = power_set(primes_low);
    const vector<mpz_class> Q_high = power_set(primes_high);

    // Inner loop temporaries
    mpz_class D, q, x_1, y_1, x_n, y_n, x_np1, y_np1, x, y;

    for (mpz_class Q_1 : Q_low) {
        // Always include p as a convince multiply into Q_1 here
        // This means q=1 is skipped but that's fine as it doesn't generate solutions.
        Q_1 *= p;

        vector<uint64_t> local_counts(omp_get_max_threads(), 0);

        #pragma omp parallel for schedule(dynamic) \
            firstprivate(local_cf_64, local_cf_128) \
            private(D, q, x_1, y_1, x_n, y_n, x_np1, y_np1, x, y)
        for (const mpz_class& Q_2 : Q_high) {
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

        }

        total[0] += Q_high.size();
        for (int i = 0; i < omp_get_max_threads(); i++) {
            total[1] += local_counts[i];
        }

        printf("\t%lu/%lu", total[1], total[0]);
    }

    if (p != P) {
        // Add all our stats to the next prime.
        p_stats[p_i + 1].combine(p_stats[p_i]);
        // Delete our found afterwards
        p_stats[p_i].found.clear();
    }
}


void verify_expand_D(char* argv1) {
    auto primes = get_primes(149);
    mpz_class D(argv1);
    vector<__uint128_t> local_cf(MAX_CF + 5, 0);
    assert(pell_solution_CF(D, local_cf));
    auto t = maybe_expand_cf(local_cf, primes);
}


int main(int argc, char** argv) {
    //verify_expand_D(argv[1]); exit(0);

    assert(argc == 2);
    assert(mp_bits_per_limb == 64);

    int n = atol(argv[1]);

    auto primes = get_primes(n);
    auto P = primes.back();
    if (n < 0 || ((unsigned) n != P)) {
        printf("Usage: %s P[=]\n", argv[0]);
        printf("P(%d) must be prime\n", n);
        exit(1);
    }

    if (1) {
        double primorial_P = std::accumulate(primes.begin(), primes.end(), 1.0, std::multiplies<>{});
        auto d = log2(primorial_P);

        // (a,b,c) = (1, n, n+1) | log(n+1) / log(rad(1*n*(n+1))) = log(n+1) / log(primorial_P)
        //
        // limit is (P^1.63)
        double limit = primorial_P + primorial_P * primorial_P;
        MAX_CF = ceil(log(limit) / log((1 + sqrt(5)) / 2));
        LIMIT_D = limit;
        LIMIT = limit;
        LIMIT_ROOT = sqrt(limit);

        if (1) { // Log of product of primes, tells us about square(2*q)
            printf("|Primes %d ... %d| = %lu -> log2(P) = %.1f | MAX_CF = %lu\n",
                   primes.front(), primes.back(), primes.size(), d, MAX_CF);
        }
    }

    vector<AllStats> p_stats;
    for (uint32_t i = 0; i < primes.size(); i++) {
        p_stats.emplace_back(primes[i], i, primes.size()-1);
    }

    bool fancy_printing = true;
    if (fancy_printing) printf("\033c"); // Clear screen

    for (uint32_t p_i = 0; p_i < primes.size(); p_i++) {
        auto p = primes[p_i];
        StormersTheorem(p, P, p_stats, fancy_printing);

        if (!fancy_printing) {
            auto& stats = p_stats[p_i];
            stats.sort_and_test_found();
            stats.print_stats(p == 2, false, p == P);
        }
    }
}
