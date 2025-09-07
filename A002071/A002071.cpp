#include <atomic>
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

#ifndef _OPENMP
    // Fakes in case -fopenmp isn't passed
    int omp_get_max_threads() { return 1; }
    int omp_get_thread_num() { return 0; }
#endif

#define XSTR(X) STR(X)
#define STR(X) #X

// 2-3x faster
#define USE_CF true

// Relates to fibonacci and max solution
size_t MAX_CF = 0;
mpz_class LIMIT = 0;
mpz_class LIMIT_ROOT = 0;

using std::atomic;
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


/*
vector<mpz_class> power_set(const vector<uint32_t>& set) {
    size_t n = set.size();
    size_t size = 1UL << n; // Equivalent to 2^n
    vector<mpz_class> ps;
    ps.reserve(size);

    for (size_t i = 1; i < size; i++) {
        mpz_class p = 1;
        for (size_t j = 0; j < n; j++) {
            if ((i >> j) & 1) {
                p *= set[j];
            }
        }
        ps.push_back(p);
    }
    return ps;
}
*/

vector<mpz_class> power_set(const vector<uint32_t>& set) {
    size_t n = set.size();
    assert( n < 30 );
    size_t size = 1 << n; // Equivalent to 2^n
    // Set all entries to 1.
    vector<mpz_class> ps;
    ps.reserve(size);
    ps.push_back(1); // For entry 0

    for (size_t i = 1; i < size; i++) {
        // Guarenteed to have some bits set
        size_t j = __builtin_ctzll(i);
        assert( i & (1 << j) );
        mpz_class p = set[j] * ps[i ^ (1 << j)];
        ps.push_back(p);
    }
    return ps;
}


bool test_smooth_small(mpz_class n, vector<uint32_t> primes) {
    if (n == 1) return true;

    mpz_class t;
    uint64_t m;

    for (auto p : primes) {
        m = mpz_fdiv_q_ui(t.get_mpz_t(), n.get_mpz_t(), p);
        while (m == 0) {
            n = t;
            if (n == 1) return true;
            m = mpz_fdiv_q_ui(t.get_mpz_t(), n.get_mpz_t(), p);
        }
    }

    return false;
}

/*
pair<mpz_class, mpz_class> pell_PQA_simple(mpz_class D) {
    // smallest solutions to x^2 - D*y^2 = 1

    mpz_class d = sqrt(D);
    mpz_class two_a0 = 2*d;

    //mpz_class P_0 = 0;
    //mpz_class Q_0 = 1;

    mpz_class A_im2 = 0;
    mpz_class A_im1 = 1;

    mpz_class B_im2 = 1;
    mpz_class B_im1 = 0;

    mpz_class G_im2 = 0;  // -P_0
    mpz_class G_im1 = 1;  // Q_0

    mpz_class P_im1, Q_im1, a_im1;

    // i = 0
    size_t i = 0;
    mpz_class P_i = 0; // P_0;
    mpz_class Q_i = 1; // Q_0;

    mpz_class a_i = d; // (P_i + d) / Q_i;
    mpz_class A_i = d; // a_i * A_im1 + A_im2;
    mpz_class B_i = 1; // a_i * B_im1 + B_im2;
    mpz_class G_i = d; // a_i * G_im1 + G_im2;

    //gmp_printf("i P Q  a A B G | %lu %Zd, %Zd   %Zd, %Zd, %Zd, %Zd\n", i, P_i, Q_i, a_i, A_i, B_i, G_i);

    // i >= 1
    for (; a_i != two_a0; ) {
        i++;
        A_im2 = A_im1;
        B_im2 = B_im1;
        G_im2 = G_im1;

        P_im1 = P_i;
        Q_im1 = Q_i;

        a_im1 = a_i;
        A_im1 = A_i;
        B_im1 = B_i;
        G_im1 = G_i;

        P_i = a_im1 * Q_im1 - P_im1;
        Q_i = (D - P_i*P_i) / Q_im1;

        a_i = (P_i + d) / Q_i;
        A_i = a_i * A_im1 + A_im2;
        B_i = a_i * B_im1 + B_im2;
        G_i = a_i * G_im1 + G_im2;
        //gmp_printf("i P Q  a A B G | %lu %Zd, %Zd   %Zd, %Zd, %Zd, %Zd\n", i, P_i, Q_i, a_i, A_i, B_i, G_i);
    }

    size_t l = i;
    assert( (a_i == two_a0) == (Q_i == 1) );

    if ((l & 1) == 0) {
        // Even Length
        return {G_im1, B_im1};
    }
    // Computer terms G_(k*l-1) from G(l-1)
    return {G_im1*G_im1 + D * B_im1*B_im1, 2 * G_im1 * B_im1};
}
*/


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
        //gmp_printf("i P Q  a A B G | %lu %Zd, %Zd   %Zd, %Zd, %Zd\n", i, P_i, Q_i, a_i, B_i, G_i);

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


inline __uint128_t from_mpz_class(const mpz_class& t) {
    // Awkward hack two limb x into t
    return (((__uint128_t) mpz_getlimbn(t.get_mpz_t(), 1)) << 64) | mpz_getlimbn(t.get_mpz_t(), 0);
}

// TODO consider using theadprivate(vector<cf> with reserved spaced)
vector<__uint128_t> continued_fraction_sqrt_128(mpz_class x_in) {
    assert( mpz_sizeinbase(x_in.get_mpz_t(), 2) <= 126 );

    mpz_class t = sqrt(x_in);
    assert( mpz_fits_ulong_p(t.get_mpz_t()) );
    __uint128_t a0 = mpz_get_ui(t.get_mpz_t());

    __uint128_t x = from_mpz_class(x_in);

    // it feels like these should work as uint64_t but they don't.
    __uint128_t b = a0;
    __uint128_t c = x - b*b;
    __uint128_t a = (a0 << 1) / c;

    // https://mathworld.wolfram.com/PeriodicContinuedFraction.html#eqn2
    // for squarefree D, 0 < ak < 2 * sqrt(n)
    std::vector<__uint128_t> cf = {a0, a};
    // TODO Check if this is faster
    //cf.reserve(MAX_CF+1);

    __uint128_t two_a0 = a0 << 1;
    for (uint32_t i = 2; i <= MAX_CF && a != two_a0; i++) {
        b = a*c - b;
        c = (x - b*b) / c;
        a = (a0 + b) / c;
        cf.push_back(a);
    }

    return cf;
}

vector<__uint128_t> continued_fraction_sqrt(mpz_class x) {
    // Assume sqrt is uint64_t
    mpz_class a0 = sqrt(x);

    mpz_class b = a0;
    mpz_class c = x - b*b;
    mpz_class a = (a0 + b) / c;

    std::vector<__uint128_t> cf = { from_mpz_class(a0), from_mpz_class(a), };

    mpz_class two_a0 = 2 * a0;
    for (uint32_t i = 2; i <= MAX_CF && a != two_a0; i++) {
        b = a*c - b;
        c = (x - b*b) / c;
        a = (a0 + b) / c;
        cf.push_back(from_mpz_class(a));
    }

    return cf;
}

vector<__uint128_t> pell_solution_CF(mpz_class n) {
    // smallest solutions to x^2 - n*y^2 = 1

    // sqrts of square free numbers are always finite.
    vector<__uint128_t> cf;
    if (mpz_sizeinbase(n.get_mpz_t(), 2) <= 126) {
        // 10x faster!
        cf = continued_fraction_sqrt_128(n);
    } else {
        cf = continued_fraction_sqrt(n);
    }
    if (cf.size() > MAX_CF)
        return {};

    assert(cf.back() == 2*cf.front());

    // https://en.wikipedia.org/wiki/Pell%27s_equation#Fundamental_solution_via_continued_fractions
    auto r = cf.size() - 1; // don't count leading value
    if (r % 2 == 0) {
        cf.pop_back();
        assert( cf.size() == r );
    } else {
        if (2*r >= MAX_CF)
            return {};
        cf.reserve(2*r);
        cf.insert(cf.end(), cf.begin() + 1, cf.end() - 1);
        assert( cf.size() == 2*r );
    }
    return cf;
}


inline pair<mpz_class, mpz_class> expand_continued_fraction_small(vector<__uint128_t>& cf) {
    assert( (cf.size() & 1) == 0 );

    mpz_class temp;
    mpz_class top = 0;
    mpz_class bottom = 1;
    for (auto v : cf | std::views::reverse) {
        top += ((uint64_t) v) * bottom;
        std::swap(top, bottom);
    }
    // Undo the last flip
    return {bottom, top};
}

inline pair<mpz_class, mpz_class> expand_continued_fraction(vector<__uint128_t>& cf) {
    if (cf[0] < (std::numeric_limits<uint64_t>::max() >> 1)) {
        return expand_continued_fraction_small(cf);
    }
    // A property of pell equation cf's
    assert( (cf.size() & 1) == 0 );

    mpz_class temp;
    mpz_class top = 0;
    mpz_class bottom = 1;
    for (auto v : cf | std::views::reverse) {
        top += ((uint64_t) v) * bottom;
        v >>= 64;
        if (v) {
            temp = ((uint64_t) v) * bottom;
            temp << 64;
            top += temp;
        }

        std::swap(top, bottom);
    }
    // Undo the last flip
    return {bottom, top};
}


inline mpz_class expand_continued_fraction_modulo32(vector<__uint128_t>& cf, uint32_t pk) {
    uint64_t top = 0;
    uint64_t bottom = 1;

    for (auto v : cf | std::views::reverse) {
        v %= pk;
        top += v * bottom;
        top %= pk;
        std::swap(top, bottom);
    }
    // Undo the last flip
    std::swap(top, bottom);
    return bottom;
}

inline mpz_class expand_continued_fraction_modulo(vector<__uint128_t>& cf, mpz_class pk) {
    mpz_class temp;
    mpz_class top = 0;
    mpz_class bottom = 1;
    for (auto v : cf | std::views::reverse) {
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


inline bool expand_continued_fraction_modulo_small(vector<__uint128_t>& cf, uint32_t p) {
    __uint128_t top = 0;
    __uint128_t bottom = 1;
    for (auto v : cf | std::views::reverse) {
        top += (v % p) * bottom;
        top %= p;
        std::swap(top, bottom);
    }
    // Undo the last flip
    std::swap(top, bottom);
    return bottom == 0;
}


inline uint32_t expand_continued_fraction_modulo_power_2(vector<__uint128_t>& cf) {
    __uint128_t top = 0;
    __uint128_t bottom = 1;
    for (auto v : cf | std::views::reverse) {
        top += (v & 0xFFFFFFFF) * bottom;
        top &= 0xFFFFFFFF;
        std::swap(top, bottom);
    }
    // Undo the last flip
    std::swap(top, bottom);
    return bottom;
}


const double PHI = (1 + sqrt(5)) / 2;
const double LOG_PHI = log2(PHI);

pair<mpz_class, mpz_class> maybe_expand_cf(vector<__uint128_t>& cf, vector<uint32_t>& primes) {
    // A continued fraction of length M is at least Fibonacci[M+1] / Fibonacci[M]
    // so y_i will be atleast ((1 + sqrt(5))/2) ^ K

    // Find highest power k, such that p^k divides the expanded continued fraction.
    // compare log2(product(p_i^k_i)) with

    if (cf.size() > MAX_CF) {
        // Is there a limit on the root solution size?
        return {-1, -1};
    }

    if (cf.size() < 60) {
        // Faster to just do it
        return expand_continued_fraction(cf);
    }
    double log_y_i = LOG_PHI * cf.size();

    double log_smooth_factors = 0;
    for (auto p : primes) {
        if (p == 2) {
            auto rem = expand_continued_fraction_modulo_power_2(cf);
            // if rem == 0, more than 32 powers of 2!
            if (rem) {
                size_t exact = 0;
                while ((rem & 1) == 0) {
                    exact += 1;
                    rem >>= 1;
                }
                // printf("\tFound %u^%lu | log2 = %.1f\n", p, exact, log(p) * exact);
                if (exact > 0) {
                    log_smooth_factors += log(p) * exact;
                }
                continue;
            }
        } else {
            // Check if p divides y_1 without GMP
            if (!expand_continued_fraction_modulo_small(cf, p)) {
                continue;
            }
        }

        // TODO might be worth handling 2 seperately as it's easy to just do 32 powers
        // and go from there.

        uint32_t k = 1;
        mpz_class p_temp = p;
        while (true) {
            k += 4;
            p_temp *= (p*p*p*p); // p < 100, p^4 < 2**27

            mpz_class m = 0;
            if (p_temp.fits_uint_p()) {
                m = expand_continued_fraction_modulo32(cf, p_temp.get_si());
            } else {
                m = expand_continued_fraction_modulo(cf, p_temp);
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
                //printf("\tFound %u^%lu | log2 = %.1f\n", p, exact, log(p) * exact);
                log_smooth_factors += log(p) * exact;
                assert(exact > 0);
                break;
            }
        }
    }

    // TODO I'd love to measure the ratio of things expanded here
    // And also the high water mark for log_smooth_factors

    // Do all the work above to prove we can skip (this print is rarely triggering at 40 terms)
    if (log_smooth_factors + 5 > log_y_i) {
        //printf("|cf| = %lu, log2(fib) > %.1f | smooth: %.1f might divide\n",
        //        cf.size(), log_y_i, log_smooth_factors);
        return expand_continued_fraction(cf);
    }

    // Can't be p-smooth
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
        AllStats(uint32_t p) : p(p) {};

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

        void print_stats(const uint64_t N, bool last) const {
            static bool header = true;
            if (header) {
                header = false;
                int l = printf("n  P"
                                      "     total max       "
                       "              total-exact max-exact "
                       "                   total1 max1      "
                       "             total1-exact max1-exact"
                       "                   total2 max2      "
                       "             total2-exact max2-exact"
                       "               total1-sqr max1-sqr  "
                       "               total1-tri max1-tri"
                       "\n");
                printf("%s\n", std::string(l - 1, '-').c_str());
            }
            gmp_printf("%-2lu %-4lu %6lu %-28Zd %6lu %-28Zd %6lu %-28Zd %6lu %-28Zd %6lu %-28Zd %6lu %-28Zd %6lu %-28Zd %6lu %Zd\n",
                N, p,
                total.count, total.max,
                total.count_exact, total.max_exact,
                total1.count, total1.max,
                total1.count_exact, total1.max_exact,
                total2.count, total2.max,
                total2.count_exact, total2.max_exact,
                total1_square.count, total1_square.max,
                total1_triangle.count, total1_triangle.max
            );
            if (last) {
                printf("\t%lu (small: %.1f%%) -> %lu (%.1f) -> %lu (%.1f) -> %lu (%.1f) -> %lu (%.1f)\n",
                        Q, 100.0 * Q_small / Q,
                        pell[0], 100.0 * pell[0] / (Q + 1e-5),
                        pell[1], 100.0 * pell[1] / (pell[0] + 1e-5),
                        pell[2], 100.0 * pell[2] / (pell[1] + 1e-5),
                        pell[3], 100.0 * pell[3] / (pell[2] + 1e-5));
            }
        }

        uint32_t p;

        // 2^primes - 1;
        uint64_t Q = 0;

        // How many of Q <= 127 Bits
        uint64_t Q_small = 0;

        // Measures early exit from MAX_CF
        // Had CF, Had possible smooth CF, had smooth fundemental, total x_n count
        uint64_t pell[4] = {};

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



AllStats StormersTheorem(vector<uint32_t> primes) {
    auto P = primes.back();
    auto solution_count = std::max<uint32_t>(3, (P + 1) / 2);

    AllStats global_count(P);
    vector<AllStats> local_counts;
    for (int i = 0; i < omp_get_max_threads(); i++) {
        local_counts.emplace_back(P);
    }

    // Minimize memory usage by breaking Q' into X chunks
    // 4 million entries * 32 bytes -> ~128 MB
    const int32_t LOW_PRIMES = primes.size() < 22 ? 0 : primes.size() - 20;
    assert( 0 <= LOW_PRIMES && (unsigned) LOW_PRIMES <= primes.size());
    assert( ((unsigned) LOW_PRIMES + 1) <= primes.size());
    vector<uint32_t> primes_low(primes.begin(), primes.begin() + LOW_PRIMES);
    vector<uint32_t> primes_high(primes.begin() + LOW_PRIMES, primes.end());
    const vector<mpz_class> Q_low = power_set(primes_low);
    const vector<mpz_class> Q_high = power_set(primes_high);
    //std::sort(Q_low.begin(), Q_low.end());
    //std::sort(Q_high.begin(), Q_high.end());

    for (mpz_class Q_1 : Q_low) {

        #pragma omp parallel for schedule(dynamic)
        for (mpz_class Q_2 : Q_high) {
            mpz_class q = Q_1 * Q_2;

            // Mucks with code, doesn't generate interesting solutions
            if (q == 1) continue;

            // Lucas computes D = q which generates other interesting numbers
            // Lehmer used D = 2 * q which only generates A002071
            mpz_class D = q; // * 2;

            AllStats &count = local_counts[omp_get_thread_num()];
            count.Q += 1;

            mpz_class x_1, y_1;
            if (USE_CF) {
                count.Q_small += (mpz_sizeinbase(D.get_mpz_t(), 2) <= 126);

                vector<__uint128_t> pell_cf = pell_solution_CF(D);
                if (pell_cf.empty()) {
                    continue;
                }

                count.pell[0] += 1;

                auto t = maybe_expand_cf(pell_cf, primes);
                x_1 = t.first;
                y_1 = t.second;

                if (y_1 < 0) {
                    // y_1 was not going to smooth
                    continue;
                }

                count.pell[1] += 1;
            } else {
                // Testing out PQa algorithm
                auto t = pell_PQA(q);
                x_1 = t.first;
                y_1 = t.second;

                if (y_1 < 0) {
                    // y_1 was not going to smooth
                    continue;
                }
                count.pell[0] += 1;
                count.pell[1] += 1;
            }

            // 1-index is better; technically solution 0 is (1, 0)
            vector<uint8_t> y_is_smooth(solution_count+1, true);

            assert( x_1 * x_1 - D * y_1 * y_1 == 1 );

            mpz_class x_n = 1;
            mpz_class y_n = 0;
            //gmp_printf("%Zd -> %Zd, %Zd\n", q, x_1, y_1);

            for (uint64_t i = 1; i <= solution_count; i++) {
                mpz_class x_np1 = x_1 * x_n + D * y_1 * y_n;
                mpz_class y_np1 = x_1 * y_n + y_1 * x_n;
                x_n = x_np1;
                y_n = y_np1;

                if (!y_is_smooth[i])
                    continue;

                count.pell[3] += 1;

                if ( mpz_even_p(D.get_mpz_t()) ) {
                    //assert( x_n * x_n - 2 * (D/2) * y_n * y_n == 1 );
                    assert( mpz_odd_p( x_n.get_mpz_t()) ); // x is always odd
                    assert( mpz_even_p(y_n.get_mpz_t()) ); // y is always even
                    // Theorem 1 (12) and (13) says y_smooth implies (x-1)/2 and (x+1)/2 are p-smooth
                }

                auto y_smooth = test_smooth_small(y_n, primes);
                if (y_smooth) {
                    if (i == 1) {
                        count.pell[2] += 1;
                    }
                    mpz_class x = x_n - 1;
                    mpz_class y = x_n + 1;

                    assert( test_smooth_small(x, primes) );
                    assert( test_smooth_small(y, primes) );

                    //gmp_printf("\t\t%Zd %Zd\n", x, y);
                    count.process_pair(x, y);
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
    }

    for (const auto& count : local_counts) {
        global_count.combine(count);
    }

    return global_count;
}

AllStats run(int n) {
    auto primes = get_primes(n);

    double P = std::accumulate(primes.begin(), primes.end(), 1.0, std::multiplies<>{});
    auto d = log2(P);

    // limit is (P + 2 * 1 * P^2)
    double limit = P + 2 * 1 * P * P;
    MAX_CF = ceil(log(limit) / log((1 + sqrt(5)) / 2));
    LIMIT = limit;
    LIMIT_ROOT = sqrt(limit);

    if (0) { // Log of product of primes, tells us about square(2*q)
        printf("Primes(%d) |%ld| log2(P) = %.1f | %d, ..., %d | MAX_CF = %lu\n",
               n, primes.size(), d, primes[0], primes.back(), MAX_CF);
    }

    return StormersTheorem(primes);
}

int main(int argc, char** argv) {
    assert(argc == 2);
    assert(mp_bits_per_limb == 64);
    //printf("USE_CF=" XSTR(USE_CF)  "\n");

    int exact = std::string(argv[1]).back() == '=';
    int n = argc <= 1 ? 47 : atol(argv[1]);
    auto primes = get_primes(n);

    uint32_t n_p = 0;
    for (auto p : primes) {
        n_p++;
        if (exact && p != primes.back()) continue;
        auto stats = run(p);
        stats.sort_and_test_found();
        stats.print_stats(n_p, p == primes.back());
    }
}

/*
Building and Verifying all CF
41 in      .2 seconds  | 869, 63927525375, 119500865855
43 in       1 seconds  | 1153, 421138799639, 574923865277
47 in       8 seconds  | 1502, 1109496723125, 2227616372734
53 in       1 minutes  | 1930, 1453579866024, 5410040985568
59 in       7 minutes  | 2454, 20628591204480, 39849760877950
61 in      31 minutes  | 3106, 31887350832896, 84714513957619
67 in 7.5/162 minutes  | 3896, 31887350832896, 107119408062061

Limiting CF to 200 -> 10^41
53 in       1 seconds  | 1930,  1453579866024,                                  5410040985568
59 in       1 seconds  | 2454,  20628591204480,                                 39849760877950
61 in       2 seconds  | 3106,  31887350832896,                                 84714513957619
67 in       3 seconds  | 3896,  31887350832896,                                 107119408062061
71 in       5 seconds  | 4839,  119089041053696,                                389016340433568
73 in       9 seconds  | 6040,  2286831727304144,                               3815472057669102
79 in      16 seconds  | 7441,  9591468737351909375,                            9597632224833515927
83 in      30 seconds  | 9179,  9591468737351909375,    17451620110781856,      9625480534208543179
89 in      56 seconds  | 11134, 9591468737351909375,    166055401586083680,     9804205815819119552
97 in       3 minutes  | 13374, 9591468737351909375,    49956990469100000,      9873978435958146488
101 in      6 minutes  | 16167, 9591468737351909375,    4108258965739505499,    14055064450964056368
103 in      12 minutes | 19507, 19316158377073923834000,19316158377073923834000,19334892305049839579209
107 in      24 minutes | 23367, 19316158377073923834000,386539843111191224
                         ^^^^^ mismatched with mathimagics, they had first 23361 then later 23372
109 in      71 minutes | 27949, 19316158377073923834000,90550606380841216610,   19565068366063659445482
113 in     143 minutes | 33233, 19316158377073923834000,205142063213188103639,  19794618108168109482594
127 in     293 minutes | 39283, 19316158377073923834000,53234795127882729824,   19895808503278115235851
131 in     604 minutes | 46166, 19316158377073923834000,4114304445616636016031, 25782826082376623386017
137 in    1273 minutes | 54150, 124225935845233319439173,124225935845233319439173,156463051279554770568122
139 in    2776 minutes | 63428, 124225935845233319439173,3482435534325338530939,165751327973248712930517
*/
