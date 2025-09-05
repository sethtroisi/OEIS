#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdint>
#include <iostream>
#include <ranges>
#include <utility>
#include <vector>

#include <gmpxx.h>


using std::vector;
using std::pair;

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


// TODO make this a generator at some point to avoid memory usage
vector<mpz_class> power_set(const vector<uint32_t>& set) {
    size_t n = set.size();
    size_t size = 1UL << n; // Equivalent to 2^n
    vector<mpz_class> ps;
    ps.reserve(size);

    for (size_t i = 0; i < size; i++) {
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


vector<uint64_t> continued_fraction_sqrt(mpz_class x) {
    // Assume sqrt is uint64_ta
    mpz_class a0 = x;
    mpz_sqrt(a0.get_mpz_t(), a0.get_mpz_t());
    assert( mpz_fits_ulong_p(a0.get_mpz_t()) );

    mpz_class b = a0;
    mpz_class c = x - b*b;
    mpz_class a = (a0 + b) / c;

    // One day I'll need to change
    assert( mpz_fits_ulong_p(a.get_mpz_t()) );
    std::vector<uint64_t> cf = {
        mpz_get_ui(a0.get_mpz_t()),
        mpz_get_ui(a.get_mpz_t()),
    };

    mpz_class two_a0 = 2 * a0;
    while (a != two_a0) {
        b = a*c - b;
        c = (x - b*b) / c;
        a = (a0 + b) / c;
        assert( mpz_fits_ulong_p(a.get_mpz_t()) );
        cf.push_back(mpz_get_ui(a.get_mpz_t()));
    }

    return cf;
}


pair<mpz_class, mpz_class> expand_continued_fraction(vector<uint64_t> cf) {
    mpz_class top = 0;
    mpz_class bottom = 1;
    for (auto v : cf | std::views::reverse) {
        top += v * bottom;
        std::swap(top, bottom);
    }
    // Undo the last flip
    return {bottom, top};
}


pair<mpz_class, mpz_class> pellSolution(mpz_class n) {
    // count smallest solutions to x^2 - n*y^2 = 1

    // sqrts of square free numbers are always finite.
    auto cf = continued_fraction_sqrt(n);

    // https://en.wikipedia.org/wiki/Pell%27s_equation#Fundamental_solution_via_continued_fractions
    auto r = cf.size() - 1; // don't count leading value
    vector<uint64_t> cf_test = cf;
    if (r % 2 == 0) {
        cf_test.pop_back();
        assert( cf_test.size() == r );
    } else {
        cf_test.reserve(2*r);
        cf_test.insert(cf_test.end(), cf.begin() + 1, cf.end() - 1);
        assert( cf_test.size() == 2*r );
    }

    return expand_continued_fraction(cf_test);
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

vector<mpz_class> StormersTheorem(vector<uint32_t> primes) {
    auto count = std::max<uint32_t>(3, (primes.back() + 1) / 2);

    vector<mpz_class> lower;
    lower.reserve(1'000'000);

    vector<mpz_class> ps = power_set(primes);
    printf("\tpowerset computed\n");

    #pragma omp parallel for schedule(dynamic, 16)
    for (mpz_class q : ps) {
        if (q == 2) {
            // Stormer doesn't like 2
            continue;
        }

        // solution 0 is (1, 0) but that's boring
        vector<uint8_t> y_is_smooth(count+1, true);

        // the count smallest solutions to x^2 - 2*q*y^2 = 1
        mpz_class n = 2 * q;
        pair<mpz_class, mpz_class> t = pellSolution(n);
        mpz_class x_1 = t.first;
        mpz_class y_1 = t.second;

        if (y_1 < 0) {
            // y_1 was not going to smooth
            continue;
        }

        assert( x_1 * x_1 - n * y_1 * y_1 == 1 );

        mpz_class x_n = 1;
        mpz_class y_n = 0;
        //gmp_printf("%Zd -> %Zd, %Zd\n", q, x_1, y_1);

        vector<mpz_class> temp;

        for (uint64_t i = 1; i <= count; i++) {

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
            	//gmp_printf("\t%Zd @ %lu -> %Zd, %Zd\n", q, i, x_n, y_n);

                assert( test_smooth_small(S, primes) );
                assert( test_smooth_small(T, primes) );
                temp.push_back(S);
            } else {
                // all future solutions y_(i*k) are divisible by y_i which is not n-smooth
                if (i == 1) {
                    // If not smooth no other solution can be smooth.
                    break;
                }
                for (size_t j = i; j <= count; j += i) {
                    y_is_smooth[j] = false;
                }
            }
        }

        if (temp.size()) {
            #pragma omp critical
            for (auto t : temp) {
                lower.push_back(t);
            }
        }
    }

    // TODO verify all unique
    return lower;
}


int main(int argc, char** argv) {
    int n = argc <= 1 ? 47 : atol(argv[1]);
    printf("\tprimes up to %d\n", n);
    auto primes = get_primes(n);
    printf("\tprimes(%d) = |%ld| = %d, %d, ..., %d\n",
           n, primes.size(), primes[0], primes[1], primes.back());


    auto M = std::max<uint32_t>(3, (primes.back() + 1) / 2);
    uint64_t maxCount = ((1L << primes.size()) - 1) * M;
    printf("\tMax possible: %lu\n", maxCount);

    auto matches = StormersTheorem(primes);
    printf("\tfound: %lu\n", matches.size());
    if (matches.size()) {
        auto [min, max] = std::minmax_element(matches.begin(), matches.end());
        gmp_printf("\tmin: %Zd\n", *min);
        gmp_printf("\tmax: %Zd\n", *max);
        printf("\n");
    }

	mpz_class summation = 0;
	for (auto m: matches) { summation += m; }
	gmp_printf("Answer 581: %Zd\n", summation);
}

/*
41 in  .4 seconds  | 869, 63927525375, 119500865855
43 in   2 seconds  | 1153, 421138799639, 574923865277
47 in  13 seconds  | 1502, 1109496723125, 2227616372734
53 in   2 minutes  | 1930, 1453579866024, 5410040985568
59 in  32 minutes  | 2454, 20628591204480, 39849760877950
*/
