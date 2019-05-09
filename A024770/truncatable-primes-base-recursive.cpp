#include <atomic>
#include <gmp.h>
#include <gmpxx.h>
#include <iostream>
#include <vector>

#pragma omp declare reduction (merge : std::vector<mpz_class> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

using namespace std;

const vector<long> SMALL_PRIMES = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37,
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

#define PRIME_REPS 25


#define LEFT false
//#define LEFT true


// Primes found with X digits
atomic<long> total;
atomic<long> depth[200] = {};

// List of primes when recurse_base stopped.
vector<mpz_class> middle_gen;

// That thing.
vector<mpz_class> left_mults;

void print_counts(int upto) {
    for (int i = 1; i <= upto & depth[i] > 0; i += (i == 1) ? 4 : 5) {
      cout << i << ":" << depth[i] << "  ";
    }
    cout << endl;
}


void recurse_base(const int base,
                  int digits,
                  bool phase_one, /* potentially helpful for gcc to see two patterns */
                  int stop_depth, /* Used to generate parallel list. */
                  mpz_class current,
                  mpz_class left_mult /* unused for LEFT = false */
) {
  if (phase_one && digits == stop_depth) {
    depth[0]++;
    middle_gen.push_back(current);
    return;
  }

  depth[digits] += 1;
  digits += 1;

#if LEFT
  left_mult *= base;
  for (int l = 1; l < base; l++) {
    mpz_class temp = l * left_mult + current;
#else

  mpz_class left = current * base;
  for (int r = 1; r < base; r += 2 - (base & 1)) {
    mpz_class temp = left + r;
#endif

    if (mpz_probab_prime_p(temp.get_mpz_t(), PRIME_REPS)) {
      total += 1;
      recurse_base(base, digits, phase_one, stop_depth, temp, left_mult);
    }
  }
}


long truncatable_primes(const int base) {
  total = 0;
  for (int i = 0; i < 200; i++) { depth[i] = 0; }
  middle_gen.clear();

  int stop_depth = 5;
  mpz_class left_mult = 1;

  // Don't syncronize as it adds to vector
  for (auto pp = SMALL_PRIMES.begin(); pp != SMALL_PRIMES.end(); pp++) {
    long p = *pp;
    if (p >= base) { continue; }

    total += 1;
    recurse_base(base, 1, true, stop_depth, p, left_mult);
  }

  if (!middle_gen.empty()) {
    cout << "\t" << total
         << " <= " << stop_depth << " leafs: " << middle_gen.size() << endl;;
  }

  for (int i = 2; i < stop_depth; i++) { left_mult *= base; }

  #pragma omp parallel for
  for (auto cur = middle_gen.begin(); cur < middle_gen.end(); cur++) {
    recurse_base(base, stop_depth, false, -1, *cur, left_mult);
  }

  cout << "\t"; print_counts(200);
  cout << endl;

  return total;
};



int
main (void)
{
  long sum = 0;
//  for (int base = 2; base <= 40; base++) {
//  for (int base = 50; base <= 53; base++) {
  for (int base = 92; base <= 100; base++) {
    long result = truncatable_primes(base);
    sum += result;
    cout << base << " " << result << endl;
  }
  cout << sum << endl;
  return 0;
}
