#include <gmp.h>
#include <gmpxx.h>
#include <iostream>
#include <vector>

#pragma omp declare reduction (merge : std::vector<mpz_class> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

using namespace std;

long SMALL_PRIMES[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37,
41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199,
211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293,
307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397,
401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491,
499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601,
607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701,
709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821,
823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929,
937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009
};

#define LEFT false
//#define LEFT true

#define PRIME_REPS 25

// TODO keep statistics about depth and print occasionally

long recurse(const int base,
  #if LEFT
             mpz_class old_left_mult,
  #endif
             mpz_class start) {
  long count = 0;

  #if LEFT
    mpz left_mult  = old_left_mult * base;
    #pragma omp parallel for reduction(+:count) if (mpz_cmp_ui(start.get_mpz_t(), 1000000) < 0) schedule(static, 1)
    for (int l = 1; l < base; l++) {
      mpz_class temp = l * left_mult + (*it_v);
      if (mpz_probab_prime_p(temp.get_mpz_t(), PRIME_REPS)) {
        count += 1 + recurse(base, left_mult, temp);
      }
    }
  #else
    mpz_class left = start * base;

    #pragma omp parallel for reduction(+:count) if (mpz_cmp_ui(start.get_mpz_t(), 1000000) < 0) schedule(static, 1)
    for (int r = 1; r < base; r += 2 - (base % 2)) {
      mpz_class temp = left + r;
      if (mpz_probab_prime_p(temp.get_mpz_t(), PRIME_REPS)) {
        count += 1 + recurse(base, temp);
      }
  #endif
  }
  return count;
}


long truncatable_primes(const int base) {
  long count = 0;

  #pragma omp parallel for reduction(+:count) schedule(static, 1)
  for (int i = 0; i < 168; i++) {
    long p = SMALL_PRIMES[i];
    if (p >= base) {
      continue;
    }

  #if LEFT
    count += 1 + recurse(base, base, p);
  #else
    count += 1 + recurse(base, p);
  #endif
  }

  return count;
};


int
main (void)
{
  for (int base = 2; base <= 60; base++) {
    long result = truncatable_primes(base);
    cout << base << " " << result << endl;
  }
  return 0;
}
