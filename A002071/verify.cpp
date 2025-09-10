#include <verify.hpp>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <vector>

#include <gmpxx.h>

using std::pair;
using std::vector;


vector<mpz_class> power_set_verify(const vector<uint32_t>& set) {
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
    std::sort(ps.begin(), ps.end());
    return ps;
}


pair<mpz_class, mpz_class> pell_PQA_verify(mpz_class D) {
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


void verify_expand_D(char* argv1) {
    auto primes = get_primes(149);
    mpz_class D(argv1);
    vector<__uint128_t> local_cf(MAX_CF + 5, 0);
    assert(pell_solution_CF(D, local_cf));
    auto t = maybe_expand_cf(local_cf, primes);
}
