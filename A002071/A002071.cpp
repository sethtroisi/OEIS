#include <atomic>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <numeric>
#include <ranges>
#include <utility>
#include <vector>

#include <gmpxx.h>


// Relates to fibonacci and max solution
#define MAX_CF 470

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


inline __uint128_t from_mpz_class(const mpz_class& t) {
    // Awkward hack two limb x into t
    return (((__uint128_t) mpz_getlimbn(t.get_mpz_t(), 1)) << 64) | mpz_getlimbn(t.get_mpz_t(), 0);
}

// TODO consider using theadprivate(vector<cf> with reserved spaced)

// For P < 101
// typedef uint64_t cf_t;
typedef __uint128_t cf_t;
vector<cf_t> continued_fraction_sqrt_128(mpz_class x_in) {
    assert( mpz_sizeinbase(x_in.get_mpz_t(), 2) < 127 );

    mpz_class t = sqrt(x_in);
    assert( mpz_fits_ulong_p(t.get_mpz_t()) );
    __uint128_t a0 = mpz_get_ui(t.get_mpz_t());

    __uint128_t x = from_mpz_class(x_in);

    __uint128_t b = a0;
    __uint128_t c = x - b*b;
    __uint128_t a = (a0 + b) / c;

    assert( a0 <= std::numeric_limits<cf_t>::max() );
    assert( a <= std::numeric_limits<cf_t>::max() );
    std::vector<cf_t> cf = {(cf_t) a0, (cf_t) a};

    __uint128_t two_a0 = 2 * a0;
    for (uint32_t i = 2; i <= MAX_CF && a != two_a0; i++) {
        b = a*c - b;
        c = (x - b*b) / c;
        a = (a0 + b) / c;
        assert( a <= std::numeric_limits<cf_t>::max() );
        cf.push_back((cf_t) a);
    }

    return cf;
}

vector<cf_t> continued_fraction_sqrt(mpz_class x) {
    // Assume sqrt is uint64_t
    mpz_class a0 = sqrt(x);

    mpz_class b = a0;
    mpz_class c = x - b*b;
    mpz_class a = (a0 + b) / c;

    std::vector<cf_t> cf = {
        (cf_t) from_mpz_class(a0),
        (cf_t) from_mpz_class(a),
    };

    mpz_class two_a0 = 2 * a0;
    for (uint32_t i = 2; i <= MAX_CF && a != two_a0; i++) {
        b = a*c - b;
        c = (x - b*b) / c;
        a = (a0 + b) / c;
        cf.push_back((cf_t) from_mpz_class(a));
    }

    return cf;
}

vector<cf_t> pell_solution_CF(mpz_class n) {
    // count smallest solutions to x^2 - n*y^2 = 1

    // sqrts of square free numbers are always finite.
    vector<cf_t> cf;
    if (mpz_sizeinbase(n.get_mpz_t(), 2) < 127) {
        // 10x faster!
        cf = continued_fraction_sqrt_128(n);
    } else {
        cf = continued_fraction_sqrt(n);
    }
    if (cf.size() >= MAX_CF)
        return {};

    // https://en.wikipedia.org/wiki/Pell%27s_equation#Fundamental_solution_via_continued_fractions
    auto r = cf.size() - 1; // don't count leading value
    if (r % 2 == 0) {
        // Technically this makes cf 1 shorter, but ignore that.
        cf.pop_back();
        assert( cf.size() == r );
    } else {
        if (2*r > MAX_CF)
            return {};
        cf.reserve(2*r);
        cf.insert(cf.end(), cf.begin() + 1, cf.end() - 1);
        assert( cf.size() == 2*r );
    }
    return cf;
}


inline mpz_class mul(uint64_t v, mpz_class &t, mpz_class &temp) {
    return v * t;
}

inline mpz_class mul(__uint128_t v, mpz_class &t, mpz_class &temp) {
   temp = t << 64;
   temp = ((uint64_t) (v >> 64)) * temp;
   temp += ((uint64_t) v) * t;
   return temp;
}


inline pair<mpz_class, mpz_class> expand_continued_fraction(vector<cf_t>& cf) {
    // A property of pell equation cf's
    assert( cf.size() % 2 == 0 );

    mpz_class temp;
    mpz_class top = 0;
    mpz_class bottom = 1;
    for (auto v : cf | std::views::reverse) {
        top += mul(v, bottom, temp);
        std::swap(top, bottom);
    }
    // Undo the last flip
    return {bottom, top};
}


inline mpz_class expand_continued_fraction_modulo(vector<cf_t>& cf, mpz_class pk) {
    mpz_class temp;
    mpz_class top = 0;
    mpz_class bottom = 1;
    for (auto v : cf | std::views::reverse) {
        top += mul(v, bottom, temp);
        top %= pk;
        std::swap(top, bottom);
    }
    // Undo the last flip
    std::swap(top, bottom);
    return bottom;
}


inline bool expand_continued_fraction_modulo_small(vector<cf_t>& cf, uint32_t p) {
    cf_t top = 0;
    cf_t bottom = 1;
    for (auto v : cf | std::views::reverse) {
        top += (v % p) * bottom;
        top %= p;
        std::swap(top, bottom);
    }
    // Undo the last flip
    std::swap(top, bottom);
    return bottom == 0;
}


inline uint32_t expand_continued_fraction_modulo_power_2(vector<cf_t>& cf) {
    cf_t top = 0;
    cf_t bottom = 1;
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

pair<mpz_class, mpz_class> maybe_expand_cf(vector<cf_t>& cf, vector<uint32_t>& primes) {
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
    //printf("|cf| = %lu, log2(fib) > %.1f\n", cf.size(), log_y_i);

    double log_smooth_factors = 0;
    for (auto p : primes) {
        if (p == 2) {
            auto rem = expand_continued_fraction_modulo_power_2(cf);
            // if rem == 0, more than 32 powers of 2!
            if (rem) {
                size_t exact = 0;
                while (rem & 1 == 0) {
                    exact += 1;
                    rem >>= 1;
                }
                // printf("\tFound %u^%lu | log2 = %.1f\n", p, exact, log(p) * exact);
                if (exact > 0) {
                    log_smooth_factors += log(p) * exact;
                }
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
        uint32_t log2_size = 0;
        mpz_class p_temp = p;
        while (true) {
            k += 4;
            p_temp *= (p*p*p*p); // p < 100, p^4 < 2**27

            //printf("\tTesting %u^%u | log2 = %u\n", p, k, log2_size);
            auto m = expand_continued_fraction_modulo(cf, p_temp);
            if (m > 0) {
                // TODO could start by removing last_log2_size count
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

vector<mpz_class> StormersTheorem(vector<uint32_t> primes) {
    auto solution_count = std::max<uint32_t>(3, (primes.back() + 1) / 2);

    atomic<uint64_t> count = 0;
    atomic<uint64_t> skips = 0;
    vector<mpz_class> found;
    found.reserve(100'000);

    // Minimize memory usage by breaking Q' into X chunks
    // 4 million entries * 32 bytes -> ~128 MB
    const uint32_t LOW_PRIMES = std::max(4, (signed) primes.size() - 22);;
    vector<uint32_t> primes_low(primes.begin(), primes.begin() + LOW_PRIMES);
    vector<uint32_t> primes_high(primes.begin() + LOW_PRIMES, primes.end());
    vector<mpz_class> Q_low = power_set(primes_low);
    for (mpz_class Q_1 : Q_low) {
        vector<mpz_class> Q_high = power_set(primes_high);
        #pragma omp parallel for schedule(dynamic)
        for (mpz_class Q_2 : Q_high) {
            mpz_class q = Q_1 * Q_2;
            if (q == 2) {
                // Stormer doesn't like 2
                continue;
            }

            // solution 0 is (1, 0) but that's boring
            vector<uint8_t> y_is_smooth(solution_count+1, true);

            // the count smallest solutions to x^2 - 2*q*y^2 = 1
            mpz_class n = 2 * q;
            vector<cf_t> pell_cf = pell_solution_CF(n);
            if (pell_cf.empty()) {
                skips += 1;
                continue;
            }

            auto t = maybe_expand_cf(pell_cf, primes);
            mpz_class x_1 = t.first;
            mpz_class y_1 = t.second;

            if (y_1 < 0) {
                // y_1 was not going to smooth
                skips += 1;
                continue;
            }

            assert( x_1 * x_1 - n * y_1 * y_1 == 1 );

            mpz_class x_n = 1;
            mpz_class y_n = 0;
            //gmp_printf("%Zd -> %Zd, %Zd\n", q, x_1, y_1);

            vector<mpz_class> temp;

            for (uint64_t i = 1; i <= solution_count; i++) {
                mpz_class x_np1 = x_1 * x_n + n * y_1 * y_n;
                mpz_class y_np1 = x_1 * y_n + y_1 * x_n;
                x_n = x_np1;
                y_n = y_np1;

                if (!y_is_smooth[i])
                    continue;

                //gmp_printf("%Zd @ %2lu -> %Zd, %Zd\n", q, i, x_n, y_n);

                assert( x_n * x_n - 2 * q * y_n * y_n == 1 );
                assert( mpz_odd_p( x_n.get_mpz_t()) ); // x is always odd
                assert( mpz_even_p(y_n.get_mpz_t()) ); // y is always even

                // Theorem 1 (12) and (13)
                auto y_smooth = test_smooth_small(y_n, primes);
                //gmp_printf("  -> %u\n", (unsigned) y_smooth);
                if (y_smooth) {
                    mpz_class S = (x_n - 1) / 2;
                    mpz_class T = S + 1;

                    assert( test_smooth_small(S, primes) );
                    assert( test_smooth_small(T, primes) );
                    temp.push_back(S);
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

            if (temp.size()) {
                #pragma omp critical
                {
                    count += temp.size();
                    for (auto t : temp) {
                        found.push_back(t);
                    }
                }
            }
        }
    }

    std::sort(found.begin(), found.end());
    auto last = std::unique(found.begin(), found.end());
    if (last != found.end()) {
        printf("\n\n\nFound Duplicates!\n\n\n");
        gmp_printf("Past the End: %Zd\n", *last);
        printf("\n\n\nRemoving %lu Duplicates!\n\n\n",
                std::distance(last, found.end()));
        found.erase(last, found.end());
    }

    printf("count: %lu\n", count.load());
    printf("skips: %lu (%.1f%%)\n", skips.load(), 100.0 * skips.load() / (1ul << primes.size()));
    return found;
}

int main(int argc, char** argv) {
    int n = argc <= 1 ? 47 : atol(argv[1]);
    auto primes = get_primes(n);
    printf("Primes(%d) = |%ld| = %d, %d, ..., %d\n",
           n, primes.size(), primes[0], primes[1], primes.back());

    { // Log of product of primes, tells us about square(2*q)
        double P = std::accumulate(primes.begin(), primes.end(), 1.0, std::multiplies<>{});
        auto d = log2(P);
        printf("\tlog2(Product(primes)) = %.1f\n", d);
    }

    { // Some debug info
        auto M = std::max<uint32_t>(3, (primes.back() + 1) / 2);
        uint64_t maxCount = ((1L << primes.size()) - 1) * M;
        mpz_class K = mpz_class::fibonacci(MAX_CF);
        double pow10 = log10(K.get_d());
        printf("Max possible count: %lu\n", maxCount);
        gmp_printf("Max possible value: K=%lu, 10^%.1f = %Zd\n", MAX_CF, pow10, K.get_mpz_t());
    }

    auto matches = StormersTheorem(primes);
    printf("\n");
    printf("found:     %lu\n", matches.size());
    if (matches.size()) {
        auto [min, max] = std::minmax_element(matches.begin(), matches.end());
        assert(*min == 1);
        gmp_printf("max:       %Zd\n", *max);
        mpz_class max_with_p = 1;
        for (auto m : matches) {
            mpz_class mp1 = m + 1;
            if (m > max_with_p && ((m % primes.back() == 0) || (mp1 % primes.back() == 0)))
                max_with_p = m;
        }
        gmp_printf("max(%lu):   %Zd\n", primes.back(), max_with_p);
    }

	mpz_class summation = 0;
	for (auto m: matches) { summation += m; }
	gmp_printf("summation: %Zd\n", summation);
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
*/
