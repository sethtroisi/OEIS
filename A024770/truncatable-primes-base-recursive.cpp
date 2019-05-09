#include <atomic>
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

// Primes found with X digits
atomic<long> total;
atomic<long> depth[200] = {};

void recurse(const int base,
             int digits,
             mpz_class current
  #if LEFT
             ,mpz_class old_left_mult
  #endif
) {

  if (total % 25000000 == 0) {
    cout << "\t" << total << "  :  ";
    for (int i = 1; i <= digits; i += (i == 1) ? 4 : 5) {
      cout << i << ":" << depth[i] << "  ";
    }
    cout << endl;
  }

  depth[digits] += 1;
  digits += 1;

  #if LEFT
    mpz left_mult  = old_left_mult * base;
    #pragma omp parallel for
    for (int l = 1; l < base; l++) {
      mpz_class temp = l * left_mult + (*it_v);
      if (mpz_probab_prime_p(temp.get_mpz_t(), PRIME_REPS)) {
        total += 1;
        recurse(base, digits, temp, left_mult);
      }
    }
  #else
    mpz_class left = current * base;

    #pragma omp parallel for
    for (int r = 1; r < base; r += 2 - (base % 2)) {
      mpz_class temp = left + r;
      if (mpz_probab_prime_p(temp.get_mpz_t(), PRIME_REPS)) {
        total += 1;
        recurse(base, digits, temp);
      }
  #endif
  }
}


long truncatable_primes(const int base) {
  total = 0;
  for (int i = 0; i < 200; i++) {
    depth[0] = 0;
  }

  for (int i = 0; i < 168; i++) {
    long p = SMALL_PRIMES[i];
    if (p >= base) {
      continue;
    }

    total += 1;

#if LEFT
    recurse(base, 1, base, p);
  #else
    recurse(base, 1, p);
  #endif
  }

  return total;
};


int
main (void)
{
  for (int base = 91; base <= 100; base++) {
    // It would be nice to have an estimate (or count of depth[5]) to generate
    // an eta.
    long result = truncatable_primes(base);
    cout << base << " " << result << endl;
  }
  return 0;
}
