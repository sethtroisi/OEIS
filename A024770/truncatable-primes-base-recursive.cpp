#include <atomic>
#include <chrono>
#include <cstdio>
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

void print_counts(int inc) {
    for (int i = 1; i <= 200 & depth[i] > 0; i += (i == 1) ? max(1, inc - 1) : inc) {
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

  auto T0 = std::chrono::high_resolution_clock::now();

  // Don't parallel, as recurse_base adds to a vector (and loop is <1s already)
  for (auto pp = SMALL_PRIMES.begin(); pp != SMALL_PRIMES.end(); pp++) {
    long p = *pp;
    if (p >= base) { continue; }

    total += 1;
    recurse_base(base, 1, true, stop_depth, p, left_mult);
  }

  auto T1 = std::chrono::high_resolution_clock::now();
  chrono::duration<double> duration = T1 - T0;

  if (!middle_gen.empty()) {
    printf("\t(%3.2fs) %7ld total, leaves(%d): %ld\n",
      duration.count(), (long)total, stop_depth, middle_gen.size());
  }

  for (int i = 2; i < stop_depth; i++) { left_mult *= base; }

  auto TLast = std::chrono::high_resolution_clock::now();

  #pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < middle_gen.size(); i++) {
    mpz_class cur = middle_gen[i];
    recurse_base(base, stop_depth, false, -1, cur, left_mult);

    auto T2 = std::chrono::high_resolution_clock::now();
    if (T2 - TLast > chrono::minutes(5)) {
      duration = T2 - T1;
      float total_minutes = duration.count() / 60.0;
      float percent = (float) i / middle_gen.size();

      // NOTE: i may not be the most recently finished middle_gen, which can
      // lead to backwards progress if print interval is to small.
      // NOTE: eta = estimated running time (an estimate of final total_minutes)
      // NOTE: eta generally overestimates during the entire run because the
      // primes are sorted in hardest to easiest order (at least for right
      // truncatable primes) because smaller leading digits are more likely to
      // be prime (e.g. 2X is more likely to be prime than 7X).
      printf("\t(%4.1fm (%d/%ld) %4.1f%%, eta %4.0fm) %10ld: ",
          total_minutes, i, middle_gen.size(), 100.0 * percent, total_minutes / percent,
          (long)total);
      print_counts(5);
      TLast = T2;
    }
  }

  cout << "\t"; print_counts(1);
  cout << endl;

  return total;
};



int
main (void)
{
  long sum = 0;
//  for (int base = 2; base <= 40; base++) {
//  for (int base = 60; base <= 63; base++) {
  for (int base = 90; base <= 100; base++) {
    long result = truncatable_primes(base);
    sum += result;
    cout << base << " " << result << endl;
  }
  cout << sum << endl;
  return 0;
}
