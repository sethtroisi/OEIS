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

#include "int128_compute_cf.h"

#ifndef _OPENMP
    // Fakes in case -fopenmp isn't passed
    int omp_get_max_threads() { return 1; }
    int omp_get_thread_num() { return 0; }
#endif

#define VERIFY true

// Relates to fibonacci and max solution
size_t MAX_CF = 10'000;
double LIMIT_D = 0;
mpz_class LIMIT = 0;
mpz_class LIMIT_ROOT = 0;

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


// Stein's Algorithm
uint64_t gcd(uint64_t a, uint64_t b)
{
    if (a == 0) return b;
    if (b == 0) return a;

    // largest power of 2 dividing both.
    int k = std::min<uint64_t>(__builtin_ctz(a), __builtin_ctz(b));
    a >>= k;
    b >>= k;

    // Remove any factors of 2 in a.
    a >>= __builtin_ctz(a);

    // a always odd below
    do
    {
        while ((b & 1) == 0)
            b >>= 1;

        // if necessary swap so that a <= b
        if (a > b)
            std::swap(a, b);

        // set b = b - a, odd - odd = even
        b = (b - a);
        b >>= 1;
    } while (b != 0);

    /* restore common factors of 2 */
    return a << k;
}

int test_smooth_small(mpz_class n, vector<uint32_t>& primes) {
    if (n == 1) { return 0; }
    assert( n > 1 ) ;

    { // Preshift out factors of 2, finding trailing zeros and right shift
        uint64_t twos = mpz_scan1(n.get_mpz_t(), 0);
        n >>= twos;
    }
    if (n == 1) { return 2; }

    mpz_class t;
    uint32_t max_p_j = 0;
    uint64_t res;
    uint64_t m = 1;
    size_t p_i = 1; // 2 already handled
    while ((p_i < primes.size()) || (m > 1)) {
        // Check if we can add another prime to m
        while (p_i < primes.size()) {
            if (__builtin_umull_overflow(m, primes[p_i], &res)) {
                break;
            } else {
                m = res;
                p_i++;
            }
        }
        if (p_i == primes.size() && m < 0xFFFF) {
            // square m giving 2 of each prime
            m = m * m;
        }

        // m is multiplication of many primes.
        res = mpz_fdiv_ui(n.get_mpz_t(), m);
        // Handles the m = 0 case
        res = m - res;
        m = gcd(res, m);
        if (m > 1) {
            //assert( 0 == mpz_fdiv_q_ui(n.get_mpz_t(), n.get_mpz_t(), m ) );
            mpz_divexact_ui(n.get_mpz_t(), n.get_mpz_t(), m);
            // Find the largest prime dividing m
            for (size_t p_j = p_i - 1; p_j > max_p_j; p_j--) {
                if (m % primes[p_j] == 0) {
                    max_p_j = p_j;
                    break;
                }
            }
            if (n == 1) {
                return primes[max_p_j];
            }
        }
    }

    return -1;
}

pair<mpz_class, mpz_class> pell_PQA(const mpz_class& D) {
    // smallest solutions to x^2 - D*y^2 = 1

    /**
     * Optimizations over pell_PQA_simple
     * 1. A_i is unused
     * 2. P_i, Q_i, a_i can be updated in that order without reference to im1
     * 3. Handle B_im2 and G_im2 with +=
     * 4. Cleaned up P_0, Q_0, B_im2, G_im2
     * 5. Have P_i, Q_i half a loop ahead
     */

    mpz_class d = sqrt(D);
    mpz_class two_a0 = 2*d;

    // Quick check that D isn't square
    if (d*d == D) return {-1, -1};

    mpz_class B_im1 = 0;
    mpz_class G_im1 = 1; // Q_0

    // i = 0
    size_t i = 0;
    mpz_class P_i = 0; // P_0;
    mpz_class Q_i = 1; // Q_0;

    mpz_class a_i = d; // (P_i + d) / Q_i;
    mpz_class B_i = 1; // a_i * B_im1 + B_im2;
    mpz_class G_i = d; // a_i * G_im1 + G_im2;

    // Advance i, P_i, Q_i and see if a, B, and G should be updated
    i++;
    P_i = a_i * Q_i - P_i;   // Q_im1, P_im1, a_im1 but here _i1
    Q_i = (D - P_i*P_i) / Q_i; // Same as above

    // i >= 1
    for (; Q_i != 1 && G_i <= LIMIT; ) {
        std::swap(B_im1, B_i);
        std::swap(G_im1, G_i);

        a_i = (P_i + d) / Q_i;
        //A_i += a_i * A_im1;
        B_i += a_i * B_im1;
        G_i += a_i * G_im1;

        i++;
        P_i = a_i * Q_i - P_i;   // Q_im1, P_im1, a_im1 but here _i1
        Q_i = (D - P_i*P_i) / Q_i; // Same as above
    }

    if (G_i > LIMIT) return {-1, -1};

    // Calc next a_i to verify
    a_i = (P_i + d) / Q_i;

    size_t l = i;
    assert( (a_i == two_a0) == (Q_i == 1) );

    // B and G are 1 index behind, which is the value needed!

    if ((l & 1) == 0) {
        // Even Length
        return {G_i, B_i};
    }

    if (G_i > LIMIT_ROOT || B_i > LIMIT_ROOT) return {-1, -1};

    // Computer terms G_(k*l-1) from G(l-1)
    return {G_i*G_i + D * B_i*B_i, 2 * G_i * B_i};
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
    __uint128_t c = x - a0*a0;
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


// https://mathworld.wolfram.com/PeriodicContinuedFraction.html#eqn2
// for squarefree D, 0 < ak < 2 * sqrt(n)
bool continued_fraction_sqrt_126(mpz_class x_in, vector<uint64_t>& cf) {
    assert( mpz_sizeinbase(x_in.get_mpz_t(), 2) <= 126 );

    mpz_class t = sqrt(x_in);
    assert( mpz_fits_ulong_p(t.get_mpz_t()) );
    __uint128_t a0 = mpz_get_ui(t.get_mpz_t());

    __uint128_t x = from_mpz_class(x_in);

    __uint128_t b = a0;
    __uint128_t c = x - b*b;
    __uint128_t a = (a0 << 1) / c;

    //cf.clear();
    //cf.push_back(a0);
    //cf.push_back(a);
    uint64_t i = 0;
    cf[++i] = a0;
    cf[++i] = a;

    __uint128_t two_a0 = a0 << 1;
    for (; i <= MAX_CF && a != two_a0; ) {
        b = a*c - b;
        c = (x - b*b) / c;
        a = (a0 + b) / c;

        // 1 <= b <= a0
        // 1 <= c <= a0 + b
        // c | (x - b*b)
        cf[++i] = a;
    }
    cf[0] = i;
    return a == two_a0;
}

bool continued_fraction_sqrt(mpz_class x, vector<__uint128_t>& cf) {
    // Assume sqrt is uint64_t
    mpz_class a0 = sqrt(x);

    mpz_class b = a0;
    mpz_class c = x - b*b;
    mpz_class a = (a0 + b) / c;

    uint64_t i = 0;
    cf[++i] = from_mpz_class(a0); // index 1
    cf[++i] = from_mpz_class(a);

    mpz_class two_a0 = 2 * a0;
    for (uint32_t i = 2; i <= MAX_CF && a != two_a0; i++) {
        b = a*c - b;
        c = (x - b*b) / c;
        a = (a0 + b) / c;
        cf[++i] = from_mpz_class(a);
    }
    cf[0] = i;
    return a == two_a0;
}

bool pell_solution_CF_126(mpz_class n, vector<uint64_t>& cf) {
    // smallest solutions to x^2 - n*y^2 = 1
    // sqrts of square free numbers are always finite.

    // 10x faster!
    assert (mpz_sizeinbase(n.get_mpz_t(), 2) <= 126);

    /*
    auto test = continued_fraction_sqrt_126_pessemistic(n);
    if (!test.first)
        return false;
    if (((test.second - 1) % 2 == 1) && (2*test.second >= MAX_CF+2))
        return false;
    */

    if (!continued_fraction_sqrt_126(n, cf)) return false;
    size_t cf_size = cf[0];
    assert( cf_size <= (MAX_CF+1) ); // Allow 1 extra so that in even case we can remove 1.
    assert( (cf[1] << 1) == cf[cf_size] );

    // https://en.wikipedia.org/wiki/Pell%27s_equation#Fundamental_solution_via_continued_fractions
    auto r = cf_size - 1; // don't count leading value
    if (r % 2 == 0) {
        cf[0]--;
    } else {
        if (2*r >= MAX_CF)
            return false;
        for (size_t i = 1; i < cf_size - 1; i++) {
            cf[cf_size + i] = cf[i+1];
        }
        cf[0] = 2*r;
    }
    return true;
}

bool pell_solution_CF(mpz_class n, vector<__uint128_t>& cf) {
    if (!continued_fraction_sqrt(n, cf)) return false;
    size_t cf_size = cf[0];
    assert( cf_size <= (MAX_CF+1) );
    assert( (cf[1] << 1) == cf[cf_size] );

    auto r = cf_size - 1;
    if (r % 2 == 0) {
        cf[0]--;
    } else {
        if (2*r >= MAX_CF)
            return false;
        for (size_t i = 1; i < cf_size - 1; i++) {
            cf[cf_size + i] = cf[i+1];
        }
        cf[0] = 2*r;
    }
    return true;
}


double expand_cf_64_as_double(vector<uint64_t>& cf) {
    size_t cf_size = cf[0];
    double temp;
    double top = 0;
    double bottom = 1;
    for (size_t i = cf_size; i > 0; i--) {
        temp = bottom * cf[i];
        top += temp;
        std::swap(top, bottom);
    }
    return bottom;
}

__attribute__((noinline))
pair<mpz_class, mpz_class> expand_cf_64(vector<uint64_t>& cf) {
    size_t cf_size = cf[0];
    assert( (cf_size & 1) == 0 );

    mpz_class temp;
    mpz_class top = 0;
    mpz_class bottom = 1;
    //for (auto v : cf | std::views::reverse) {
    for (size_t i = cf_size; i > 0; i--) {
        auto v = cf[i];
        temp = bottom;
        temp *= v;
        top += temp;
        std::swap(top, bottom);
    }
    // Undo the last flip
    return {bottom, top};
}

__attribute__((noinline))
pair<mpz_class, mpz_class> expand_cf_small(vector<__uint128_t>& cf) {
    size_t cf_size = cf[0];
    assert( (cf.size() & 1) == 0 );

    mpz_class temp;
    mpz_class top = 0;
    mpz_class bottom = 1;
    for (size_t i = cf_size; i > 0; i--) {
        uint64_t v = cf[i];
        temp = bottom;
        temp *= v;
        top += temp;
        std::swap(top, bottom);
    }
    // Undo the last flip
    return {bottom, top};
}

pair<mpz_class, mpz_class> expand_cf(vector<__uint128_t>& cf) {
    if (cf[0] < (std::numeric_limits<uint64_t>::max() >> 2)) {
        assert(false && "Should not happen very often anymore");
        return expand_cf_small(cf);
    }
    // A property of pell equation cf's
    size_t cf_size = cf[0];
    assert( (cf_size & 1) == 0 );

    mpz_class temp;
    mpz_class top = 0;
    mpz_class bottom = 1;
    for (size_t i = cf_size; i > 0; i--) {
        auto v = cf[i];
        temp = bottom;
        temp *= (uint64_t) v;
        top += temp;
        v >>= 64;
        if (v) {
            temp = bottom;
            temp *= (uint64_t) v;
            temp << 64;
            top += temp;
        }

        std::swap(top, bottom);
    }
    // Undo the last flip
    return {bottom, top};
}


__attribute__((noinline))
uint32_t expand_cf_64_modulo32(vector<uint64_t>& cf, uint32_t pk) {
    uint64_t top = 0;
    uint64_t bottom = 1;

    size_t cf_size = cf[0];
    for (size_t i = cf_size; i > 0; i--) {
        auto v = cf[i];
        uint32_t v_mod = (v < pk) ? v : (v % pk);
        top += v_mod * bottom;
        top %= pk;
        std::swap(top, bottom);
    }
    // Undo the last flip
    std::swap(top, bottom);
    return bottom;
}

__attribute__((noinline))
uint32_t expand_cf_modulo32(vector<__uint128_t>& cf, uint32_t pk) {
    uint64_t top = 0;
    uint64_t bottom = 1;

    size_t cf_size = cf[0];
    for (size_t i = cf_size; i > 0; i--) {
        auto v = cf[i];
        uint32_t v_mod = (v < pk) ? v : (v % pk);
        top += v_mod * bottom;
        top %= pk;
        std::swap(top, bottom);
    }
    // Undo the last flip
    std::swap(top, bottom);
    return bottom;
}


__attribute__((noinline))
uint64_t expand_cf_64_modulo64(vector<uint64_t>& cf, uint64_t pk) {
    uint64_t top = 0;
    uint64_t bottom = 1;

    size_t cf_size = cf[0];
    for (size_t i = cf_size; i > 0; i--) {
        auto v = cf[i];
        if (v > pk) {
            v %= pk;
        }
        __uint128_t temp = bottom;
        temp *= v;
        temp += top;
        temp %= pk;

        top = bottom;
        bottom = (uint64_t) temp;
    }
    // Undo the last flip
    std::swap(top, bottom);
    return bottom;
}

__attribute__((noinline))
uint64_t expand_cf_modulo64(vector<__uint128_t>& cf, uint64_t pk) {
    uint64_t top = 0;
    uint64_t bottom = 1;

    size_t cf_size = cf[0];
    for (size_t i = cf_size; i > 0; i--) {
        auto v = cf[i];
        // ((v % pk) * bottom) % pk
        if (v > pk) {
            v %= pk;
        }
        // bottom * v is uint64 * uint64 might be possible to improve over uint128*uint128
        __uint128_t temp = bottom;
        temp *= (uint64_t) v;
        temp += top;
        temp %= pk;

        top = bottom;
        bottom = (uint64_t) temp;
    }
    // Undo the last flip
    std::swap(top, bottom);
    return bottom;
}


mpz_class expand_cf_modulo(vector<__uint128_t>& cf, mpz_class pk) {
    mpz_class temp;
    mpz_class top = 0;
    mpz_class bottom = 1;
    size_t cf_size = cf[0];
    for (size_t i = cf_size; i > 0; i--) {
        auto v = cf[i];
        top += ((uint64_t) v) * bottom;
        v >>= 64;
        if (v) {
            temp = ((uint64_t) v) * bottom;
            temp << 64;
            top += temp;
        }
        top %= pk;
        std::swap(top, bottom);
    }
    // Undo the last flip
    std::swap(top, bottom);
    return bottom;
}

bool expand_cf_modulo_small(vector<__uint128_t>& cf, uint32_t p) {
    __uint128_t top = 0;
    __uint128_t bottom = 1;
    size_t cf_size = cf[0];
    for (size_t i = cf_size; i > 0; i--) {
        auto v = cf[i];
        if (v < p) {
            uint64_t temp = v;
            temp *= bottom;
            temp += (uint64_t) top; // can't overflow
            temp %= p; // uint64_t is much faster
            top = temp;
        } else {
            top += (v % p) * bottom;
            top %= p;
        }
        std::swap(top, bottom);
    }
    // Undo the last flip
    std::swap(top, bottom);
    return bottom == 0;
}


/** compute expand(cf) % 2^64 */
uint32_t expand_cf_modulo_power_2(vector<__uint128_t>& cf) {
    uint64_t top = 0;
    uint64_t bottom = 1;
    size_t cf_size = cf[0];
    for (size_t i = cf_size; i > 0; i--) {
        auto v = cf[i];
        top += ((uint64_t) v) * bottom;
        std::swap(top, bottom);
    }
    // Undo the last flip
    std::swap(top, bottom);
    return bottom;
}

/** compute expand(cf) % 2^64 */
uint32_t expand_cf_64_modulo_power_2(vector<uint64_t>& cf) {
    uint64_t top = 0;
    uint64_t bottom = 1;
    size_t cf_size = cf[0];
    for (size_t i = cf_size; i > 0; i--) {
        auto v = cf[i];
        top += v * bottom;
        std::swap(top, bottom);
    }
    // Undo the last flip
    std::swap(top, bottom);
    return bottom;
}


double compute_smooth_size_verify(vector<__uint128_t>& cf, vector<uint32_t>& primes) {
    double log_smooth_factors = 0;
    for (auto p : primes) {
        if (p == 2) {
            auto rem = expand_cf_modulo_power_2(cf);
            // if rem == 0, more than 32 powers of 2!
            if (rem) {
                size_t exact = 0;
                while ((rem & 1) == 0) {
                    exact += 1;
                    rem >>= 1;
                }
                if (exact > 0) {
                    log_smooth_factors += log(p) * exact;
                    //printf("\tFound1  %u^%lu\n", p, exact);
                }
                continue;
            }
        } else {
            // Check if p divides y_1 without GMP
            if (!expand_cf_modulo_small(cf, p)) {
                continue;
            }
        }
        uint32_t k = 1;
        mpz_class p_temp = p;
        while (true) {
            k += 4;
            p_temp *= (p*p*p*p); // p < 100, p^4 < 2**27

            mpz_class m = 0;
            if (p_temp.fits_uint_p()) {
                m = expand_cf_modulo32(cf, p_temp.get_si());
            } else {
                m = expand_cf_modulo(cf, p_temp);
            }
            if (m > 0) {
                // TODO could start by removing previous k count
                size_t exact = 0;
                auto r = mpz_fdiv_q_ui(m.get_mpz_t(), m.get_mpz_t(), p);
                while (r == 0) {
                    exact += 1;
                    if (m < p) break;
                    r = mpz_fdiv_q_ui(m.get_mpz_t(), m.get_mpz_t(), p);
                }
                //printf("\tFound1  %u^%lu\n", p, exact);
                log_smooth_factors += log(p) * exact;
                assert(exact > 0);
                break;
            }
        }
    }
    return log_smooth_factors;
}

/** Handle the harder case where p is needed to large power */
uint32_t count_prime_power_in_expanded(bool is_small, vector<uint64_t>& cf_64, vector<__uint128_t>& cf, uint32_t prime) {
    // Can get 4 powers without checking for overflow
    assert( prime <= 255 );

    // largest power of prime that fits in uint32
    uint64_t power = prime * prime;
    power = power * power;

    uint64_t t = power * prime;
    while (t < std::numeric_limits<uint32_t>::max()) {
        power = t;
        t *= prime;
    }

    uint64_t m = is_small ? expand_cf_64_modulo32(cf_64, power) : expand_cf_modulo32(cf, power);
    if (m == 0) {
        // Might be able to fit one more prime in uint64_t.
        power = power * power;
        // This never overflowed at p=151 so no need yet to add that logic.
        m = is_small ? expand_cf_64_modulo64(cf_64, power) : expand_cf_modulo64(cf, power);
    }
    assert( m != 0 && "Assume never happens, can implement mpz_class if needed" );

    // count of prime that divides m
    uint32_t k = 0;
    while (true) {
        uint64_t d = m / prime;
        uint64_t r = m - d * prime;
        if (r != 0) break;
        k += 1;
        m = d;
    }
    return k;
}

double compute_smooth_size(bool is_small, vector<uint64_t>& cf_64, vector<__uint128_t>& cf, vector<uint32_t>& primes) {
    double log_smooth_factors = 0;

    { // Handle 2 first
        auto rem = is_small ? expand_cf_64_modulo_power_2(cf_64) : expand_cf_modulo_power_2(cf);
        if (rem) {
            uint32_t k = __builtin_ctz(rem);
            assert( rem & (1 << k) );
            if (k > 0) {
                log_smooth_factors += log(2) * k;
                //printf("\tFound2  %u^%u\n", 2, k);
            }
        } else {
            // if rem == 0, more than 64 powers of 2!
            assert(false && "Assume never happens, can implement if needed");
        }
    }

    size_t p_i = 1; // 2 already handled

    const vector<pair<uint32_t, int32_t>> groups_many = {
        {3*3*3*3*3*3*3 * 5*5*5*5*5 * 7*7*7u, 3},
        {11*11 * 13*13 * 17*17 * 19*19u, 4},
        {23*23 * 29*29 * 31*31u, 3},
    };

    // Handle small primes to powers.
    for (auto [K, K_count] : groups_many) {
        K_count = std::min<int32_t>(K_count, primes.size() - p_i);
        auto m_combined = is_small ? expand_cf_64_modulo32(cf_64, K) : expand_cf_modulo32(cf, K);
        for (; K_count > 0; K_count--, p_i++ ) {
            uint32_t prime = primes[p_i];
            uint32_t m = m_combined % prime;
            if (m == 0) {
                uint32_t k = 2;
                uint32_t p_k = prime * prime;
                assert (K % p_k == 0 ); // Can clean up or disable at some point
                do {
                    m = m_combined % p_k;
                    if (m != 0) {
                        k -= 1;
                        break;
                    }
                    k++;
                    p_k *= prime;
                } while (K % p_k == 0);

                if (m == 0) {
                    // around 1% of the time.
                    k = count_prime_power_in_expanded(is_small, cf_64, cf, prime);
                }
                if (k == 0) {
                    auto t = is_small ? expand_cf_64(cf_64) : expand_cf(cf);
                    gmp_printf("ERROR WITH FACTORING: %Zd, %u part of %u, %u\n", t.second, prime, K, m_combined);
                }
                assert( k > 0 );
                log_smooth_factors += log(prime) * k;
                //printf("\tFound2  %u^%u\n", prime, k);
            }
        }
        if (p_i == primes.size())
            return log_smooth_factors;
    }

    // multiplication, number of primes
    // TODO could probably reduce by 1 around 109 by rearranging order of factors. Would require more code.
    const vector<pair<uint32_t, int32_t>> groups = {
        {37*41*43*47*53u, 5},   // 27.3 bits
        {59*61*67*71*73u, 5},   // 30.2 bits
        {79*83*89*97u, 4},      // 25.8 bits
        {101*103*107*109u, 4},  // 26.9 bits
        {113*127*131*137u, 4},  // 27.9 bits
        {139*149*151*157u, 4},  // 28.9 bits
        {163*167*173*179u, 4},  // 29.7 bits
        {181*191*193*197u, 4},  // 30.3 bits
    };

    for (auto [K, K_count] : groups) {
        auto m_combined = is_small ? expand_cf_64_modulo32(cf_64, K) : expand_cf_modulo32(cf, K);
        K_count = std::min<int32_t>(K_count, primes.size() - p_i);
        for (; K_count > 0; K_count--, p_i++ ) {
            uint32_t prime = primes[p_i];
            uint32_t m = m_combined % prime;
            if (m == 0) {
                // prime appears to k power in expanded cf.
                auto k = count_prime_power_in_expanded(is_small, cf_64, cf, prime);
                if (k == 0) {
                    auto t = is_small ? expand_cf_64(cf_64) : expand_cf(cf);
                    gmp_printf("ERROR WITH FACTORING: %Zd, %u part of %u, %u\n", t.second, prime, K, m_combined);
                }
                assert( k > 0 );
                log_smooth_factors += log(prime) * k;
                //printf("\tFound2  %u^%u\n", prime, k);
            }
        }
        if (p_i == primes.size())
            return log_smooth_factors;
    }
    assert(p_i == primes.size());
    return log_smooth_factors;
}

const double PHI = (1 + sqrt(5)) / 2;
const double LOG_PHI = log2(PHI);

pair<mpz_class, mpz_class> maybe_expand_cf_64(vector<uint64_t>& cf, vector<__uint128_t>& temp, vector<uint32_t>& primes) {
    /**
     * A continued fraction of length M is at least Fibonacci[M+1] / Fibonacci[M]
     * so y_i will be atleast ((1 + sqrt(5))/2) ^ K
     *
     * Find highest power k, such that p^k divides the expanded continued fraction.
     * compare log2(product(p_i^k_i)) with fibonacci(|cf|)
     */

    size_t cf_size = cf[0];
    if (cf_size > MAX_CF) {
        return {-1, -1};
    }

    // TODO could try and expand with doubles to see if > LIMIT_D

    if (cf_size < 55) {
        return expand_cf_64(cf);
    }

    double log_y_i = LOG_PHI * cf_size;

    double log_smooth_factors = compute_smooth_size(true, cf, temp, primes);
    //double log_smooth_factors_alt = compute_smooth_size(false, cf, temp, primes);
    //assert( log_smooth_factors == log_smooth_factors_alt );
    //log_smooth_factors_alt = compute_smooth_size_verify(temp, primes);
    //assert( abs(log_smooth_factors - log_smooth_factors_alt) < 1e-4 );

    if (log_smooth_factors + 1 > log_y_i) {
        return expand_cf_64(cf);
    }

    // Can't be p-smooth
    return {-1, -1};
}

pair<mpz_class, mpz_class> maybe_expand_cf(vector<__uint128_t>& cf, vector<uint32_t>& primes) {
    size_t cf_size = cf[0];
    if (cf_size > MAX_CF) {
        return {-1, -1};
    }

    if (cf_size < 50) {
        return expand_cf(cf);
    }
    double log_y_i = LOG_PHI * cf_size;

    vector<uint64_t> temp;
    double log_smooth_factors = compute_smooth_size(false, temp, cf, primes);
    //double log_smooth_factors_alt = compute_smooth_size_verify(cf, primes);
    //assert( abs(log_smooth_factors - log_smooth_factors_alt) < 1e-4 );

    if (log_smooth_factors + 1 > log_y_i) {
        return expand_cf(cf);
    }
    return {-1, -1};
}


class StatCount {
    public:
        uint64_t  count = 0;
        mpz_class max = 0;

        uint64_t  count_exact = 0;
        mpz_class max_exact = 0;

        void process(mpz_class x, bool is_exact) {
            count += 1;
            if (x > max)
                max = x;

            if (is_exact) {
                count_exact += 1;
                if (x > max_exact) {
                    max_exact = x;
                }
            }
        }

        void combine(const StatCount& other, bool same_exact) {
            count += other.count;
            if (other.max > max) {
                max = other.max;
            }

            if (same_exact) {
                count_exact += other.count_exact;
                if (other.max_exact > max_exact) {
                    max_exact = other.max_exact;
                }
            }
        }
};

class AllStats {
    public:
        AllStats(uint32_t p, uint32_t p_i, uint32_t goal_i) : p(p), p_i(p_i), goal_i(goal_i) {};

        void process_pair(mpz_class x, mpz_class y) {
            found.push_back(x);

            uint32_t x_mod_p = mpz_fdiv_ui(x.get_mpz_t(), p);

            bool p_exact = (x_mod_p == 0) || (x_mod_p == (p - 2));
            total.process(x, p_exact);

            bool x_even = mpz_even_p( x.get_mpz_t() );
            bool y_even = mpz_even_p( y.get_mpz_t() );
            if (x_even && y_even) {
                total1.process(x/2, p_exact);


                mpz_class y_root = sqrt(y/2);
                if (2 * y_root * y_root == y) {
                    total1_square.process(x/2, p_exact);
                }
                y_root = sqrt(y);
                if (y_root * (y_root + 1) == y) {
                    total1_triangle.process(x/2, p_exact);
                }

            } else {
                total2.process(x, p_exact);
            }
        }

        void combine(const AllStats &other) {
            bool same_exact = (p == other.p);
            total.combine(other.total, same_exact);
            total1.combine(other.total1, same_exact);
            total2.combine(other.total2, same_exact);
            total1_square.combine(other.total1_square, same_exact);
            total1_triangle.combine(other.total1_triangle, same_exact);

            found.insert(found.end(), other.found.begin(), other.found.end());

            Q += other.Q;
            Q_small += other.Q_small;
            pell[0] += other.pell[0];
            pell[1] += other.pell[1];
            pell[2] += other.pell[2];
            pell[3] += other.pell[3];

            if (longest_cf < other.longest_cf) {
                longest_cf = other.longest_cf;
                longest_D = other.longest_D;
            }
        }

        void sort_and_test_found() {
            std::sort(found.begin(), found.end());
            auto last = std::unique(found.begin(), found.end());
            if (last != found.end()) {
                printf("\n\n\nFound Duplicates!\n\n\n");
                gmp_printf("Past the End: %Zd\n", *last);
                printf("\n\n\nRemoving %lu Duplicates!\n\n\n",
                        std::distance(last, found.end()));
                found.erase(last, found.end());
            }
        }

        void print_stats(bool header, bool add_dashes, bool last) const {
            if (header) {
                printf("A002071 / A002072 Solver by Seth Troisi\n");
                int l = printf("n  P"
                                    "     total max       "
                       "            total-exact max-exact "
                       "                 total1 max1      "
                       "           total1-exact max1-exact"
                       "                 total2 max2      "
                       "           total2-exact max2-exact"
                       "             total1-sqr max1-sqr  "
                       "             total1-tri max1-tri  "
             //            "              longest CF D         "
             //          "                 largest X_0           "
                       "\n");
                printf("%s\n", std::string(l - 1, '-').c_str());
            }
            auto width = gmp_printf("%-2lu %-4u %6lu %-26Zd %6lu %-26Zd %6lu %-26Zd %6lu %-26Zd %6lu %-26Zd %6lu %-26Zd %6lu %-26Zd %6lu %-26Zd\n",
            //auto width = gmp_printf("%-2lu %-4u %6lu %-26Zd %6lu %-26Zd %6lu %-26Zd %6lu %-26Zd %6lu %-26Zd %6lu %-26Zd %6lu %-26Zd %6lu %-26Zd %8lu %-26Zd\n",
                p_i, p,
                total.count, total.max,
                total.count_exact, total.max_exact,
                total1.count, total1.max,
                total1.count_exact, total1.max_exact,
                total2.count, total2.max,
                total2.count_exact, total2.max_exact,
                total1_square.count, total1_square.max,
                total1_triangle.count, total1_triangle.max
                //,longest_cf, longest_D
            );
            if (add_dashes) {
                // would be nice for this to be above the line it's related to, but width isn't exactly know yet.
                uint64_t total_Q = (2LL << p_i) - 1;
                // Happens if we primes were skipped for some reason
                uint64_t finished_Q = 2*Q > total_Q ? 2*Q - total_Q : Q;
                char buffer[1000];
                auto text_size = snprintf(
                        buffer, sizeof(buffer),
                        "  Working on P: %-2u ---- Status: %lu/%lu (%.1f%%, %.2f%%)  ",
                        p, Q, total_Q,
                        100.0 * finished_Q / total_Q, 100.0 * Q / (2LL << goal_i));
                int dashes = std::max<int>(0, (width - 1) - text_size);
                int shift = std::max<int>(1, finished_Q * dashes / total_Q);
                printf("%s%s%s\n",
                        std::string(shift, '-').c_str(),
                        buffer,
                        std::string(dashes - shift, '-').c_str());
            }
            if (last) {
                printf("\t%lu (small: %.1f%%) -> %lu (%.1f) -> %lu (%.1f) -> %lu (%.1f) -> %lu (%.1f) | MAX_CF: %lu\n",
                        Q, 100.0 * Q_small / Q,
                        pell[0], 100.0 * pell[0] / (Q + 1e-5),
                        pell[1], 100.0 * pell[1] / (pell[0] + 1e-5),
                        pell[2], 100.0 * pell[2] / (pell[1] + 1e-5),
                        pell[3], 100.0 * pell[3] / (pell[2] + 1e-5),
                        MAX_CF);
            }
        }

        // prime and prime_index and the index of the goal prime.
        uint32_t p;
        uint32_t p_i;
        uint32_t goal_i;

        // 2^primes - 1;
        uint64_t Q = 0;

        // How many of Q <= 127 Bits
        uint64_t Q_small = 0;

        // Measures early exit from MAX_CF
        // Had CF, Had possible smooth CF, had smooth fundemental, total x_n count
        uint64_t pell[4] = {};

        // Pell Equation max stats.
        uint64_t longest_cf = 0;
        mpz_class longest_D = 0;

        // x, x+2 are p-smooth
        StatCount total;

        // x, x+1 are p-smooth (AKA x, x+2 are even)
        StatCount total1;

        // x, x+2 are p-smooth (and x or x+2 is odd) (AKA not from (x, x+1) pair)
        StatCount total2;

        // x, x+1 are p-smooth and x+1 is square
        StatCount total1_square;

        // x, x+1 are p-smooth and x+1 is triangule number
        StatCount total1_triangle;


        vector<mpz_class> found;
};


/**
 * Given prime p, P with p <= P
 * handle Q = powerset([2, 3, 7, 11, p])
 */
void StormersTheorem(uint32_t p, uint32_t P, vector<AllStats>& p_stats, bool fancy_printing) {
    auto primes = get_primes(P);
    assert( P == primes.back() );
    auto solution_count = std::max<uint32_t>(3, (P + 1) / 2);

    vector<int> p_index(P+1, 0);
    for (size_t i = 0; i < primes.size(); i++) p_index[primes[i]] = i;

    int p_i = p_index[p];
    assert( primes[p_i] == p );

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

    // GPU variables.
    vector<pair<__uint128_t, __uint128_t>> temp_Q(Q_high.size(), {0, 0});
    vector<uint32_t> gpu_valid(Q_high.size(), 0);
    PessemisticCf gpu_tester(Q_high.size());

    uint64_t local_max_cf = std::min<uint64_t>(MAX_CF + 5, 1ULL << (2 * p_i + 1));
    // Inner loop temporaries
    mpz_class D, q, x_1, y_1, x_n, y_n, x_np1, y_np1, x, y, t;
    vector<uint64_t> local_cf_64(local_max_cf, 0);
    vector<__uint128_t> local_cf_128(local_max_cf, 0);

    for (mpz_class Q_1 : Q_low) {
        // Always include p as a convince multiply into Q_1 here
        // This means q=1 is skipped but that's fine as it doesn't generate solutions.
        Q_1 *= p;

        vector<vector<AllStats>> local_counts(omp_get_max_threads());
        for (int i = 0; i < omp_get_max_threads(); i++) {
            for (auto s : p_stats) {
                local_counts[i].emplace_back(s.p, s.p_i, s.goal_i);
            }
        }

        #pragma omp parallel for schedule(dynamic) private(D, t)
        for (size_t i = 0; i < Q_high.size(); i++) {
            const mpz_class& Q_2 = Q_high[i];
            D = Q_1 * Q_2;
            t = sqrt(D);
            temp_Q[i] = {from_mpz_class(D), from_mpz_class(t)};
        }

        gpu_tester.run(MAX_CF, temp_Q, gpu_valid, false);

        #pragma omp parallel for schedule(dynamic) \
            firstprivate(local_cf_64, local_cf_128) \
            private(D, q, x_1, y_1, x_n, y_n, x_np1, y_np1, x, y)
        //for (const mpz_class& Q_2 : Q_high) {
        for (size_t i = 0; i < Q_high.size(); i++) {
            q = Q_1 * Q_high[i];
            //q = Q_1 * Q_2;

            // Lucas computes D = q which generates other interesting numbers
            // Lehmer used D = 2 * q which only generates A002071
            D = q; // * 2;

            AllStats &count = local_counts[omp_get_thread_num()][p_i];
            count.Q += 1;
            bool is_small = (mpz_sizeinbase(D.get_mpz_t(), 2) <= 126);
            count.Q_small += is_small;

            if (!gpu_valid[i]) continue;

            bool valid = is_small
              ? pell_solution_CF_126(D, local_cf_64)
              : pell_solution_CF(D, local_cf_128);

            size_t cf_size = is_small ? local_cf_64[0] : local_cf_128[0];
            if (cf_size > count.longest_cf) {
                 // This can be "doubled" the value from the CF[sqrt[D]]
                 count.longest_cf = cf_size;
                 count.longest_D = D;
            }

            if (!valid) continue;

            count.pell[0] += 1;

            auto t = is_small
                ? maybe_expand_cf_64(local_cf_64, local_cf_128, primes)
                : maybe_expand_cf(local_cf_128, primes);
            x_1 = t.first;
            y_1 = t.second;

            if (y_1 < 0) {
                // y_1 was not going to smooth
                continue;
            }

            count.pell[1] += 1;

            if (x_1 < LIMIT) {
                auto t = pell_PQA(q);
                if ( !((x_1 == t.first) && (y_1 == t.second)) ) {
                    gmp_printf("Mismatch solving Pell %Zd -> (%Zd, %Zd) vs (%Zd, %Zd)\n", x_1, y_1, t.first, t.second);
                }
            }

            // 1-index is better; technically solution 0 is (1, 0)
            vector<uint8_t> y_is_smooth(solution_count+1, true);

            assert( x_1 * x_1 - D * y_1 * y_1 == 1 );

            x_n = 1;
            y_n = 0;
            //gmp_printf("Pell solution: %Zd -> %Zd, %Zd\n", D, x_1, y_1);

            for (uint64_t i = 1; i <= solution_count; i++) {
                x_np1 = x_1 * x_n + D * y_1 * y_n;
                y_np1 = x_1 * y_n + y_1 * x_n;
                x_n = x_np1;
                y_n = y_np1;

                if (!y_is_smooth[i])
                    continue;

                count.pell[3] += 1;

                //gmp_printf("\tPell solution(%lu): %Zd -> %Zd, %Zd\n", i, D, x_n, y_n);
                if ( mpz_even_p(D.get_mpz_t()) ) {
                    // Theorem 1 (12) and (13) says y_smooth implies (x-1)/2 and (x+1)/2 are p-smooth
                    //assert( x_n * x_n - 2 * (D/2) * y_n * y_n == 1 );
                    assert( mpz_odd_p( x_n.get_mpz_t()) ); // x is always odd
                    assert( mpz_even_p(y_n.get_mpz_t()) ); // y is always even
                }

                /* p-smooth(y_n) or -1 if y_n is not P-smooth */
                auto y_smooth = test_smooth_small(y_n, primes);
                /*auto y_smooth_verify = test_smooth_small_verify(y_n, primes);
                if (y_smooth != y_smooth_verify) {
                    gmp_printf("smooth not matching: %Zd\n", y_n);
                    assert( false );
                }*/

                if (y_smooth >= 0) {
                    if (i == 1) {
                        //gmp_printf("Pell solution(%lu): %Zd -> %Zd, %Zd\n", i, D, x_n, y_n);
                        count.pell[2] += 1;
                    }
                    x = x_n - 1;
                    y = x_n + 1;

                    auto a_smooth = test_smooth_small(x, primes);
                    auto b_smooth = test_smooth_small(y, primes);
                    assert( a_smooth <= std::max<int>(p, y_smooth) );
                    assert( b_smooth <= std::max<int>(p, y_smooth) );

                    auto min_smooth = std::max(a_smooth, b_smooth);
                    //gmp_printf("Result %Zd %lu -> %Zd %Zd | %d-smooth -> %d %d -> index: %d\n",
                    //        D, i, x_n, y_n, y_smooth, a_smooth, b_smooth, p_index[min_smooth]);

                    // Process to the correct P for this (x, y)
                    auto& stat = local_counts[omp_get_thread_num()][p_index[min_smooth]];
                    stat.process_pair(x, y);
                } else {
                    // all future solutions y_(i*k) are divisible by y_i which is not n-smooth
                    if (i == 1) {
                        // If not smooth no other solution can be smooth.
                        break;
                    }
                    for (size_t j = i; j <= solution_count; j += i) {
                        y_is_smooth[j] = false;
                    }
                }
            }
        }

        // Combine all the local stats for each p
        for (size_t p_j = 0; p_j < primes.size(); p_j++) {
            for (int i = 0; i < omp_get_max_threads(); i++) {
                p_stats[p_j].combine(local_counts[i][p_j]);
            }
        }

        if (fancy_printing) {
            auto lastQ = Q_1 == p * Q_low.back();
            auto last = lastQ && (p == P);
            printf("\033[H"); // Go to home position
            for (size_t p_i = 0; p_i < primes.size(); p_i++) {
                auto t = primes[p_i];
                p_stats[p_i].print_stats(t == 2, t == p, last && t == P);
            }
            // Last isn't rolled forward so manually print
            if (!last) {
                const auto& s = p_stats[p_i];
                printf("\t%lu (small: %.1f%%) -> %lu (%.1f)\n",
                        s.Q, 100.0 * s.Q_small / s.Q,
                        s.pell[0], 100.0 * s.pell[0] / (s.Q + 1e-5));
            }
        }
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
