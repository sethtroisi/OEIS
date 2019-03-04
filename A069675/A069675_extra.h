// Note the omp parallel for requires -fopenmp at compile time
// yields ~10x speedup.

#include <gmpxx.h>

#include <cassert>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>

#include "A069675_sieve.h"


using namespace std;

void AssertDivisible(int a, int d, int b, long p) {
  mpz_class t = 10;
  mpz_class mod = p;
  mpz_powm_ui(t.get_mpz_t(), t.get_mpz_t(), d, mod.get_mpz_t());

  t *= a;
  t += b;

  t %= mod;

  if (t != 0) {
    cout << "WHAT: " << a << " " << d << " " << b << " " << p << endl;
  }
  assert(t == 0);
}

void FilterSimple() {
  // Filter divisible by 2 and 5 mods
  for (int a = 1; a <= 9; a += 1) {
    for (int test_d = 1; test_d <= MAX_DIGITS; test_d += 1) {
      is_prime[test_d][a][2] = 2;
      is_prime[test_d][a][4] = 2;
      is_prime[test_d][a][6] = 2;
      is_prime[test_d][a][8] = 2;
      is_prime[test_d][a][5] = 5;
    }
  }

  // Filter divisible by 3 mods.
  for (int a = 1; a <= 9; a += 1) {
    for (int b = 1; b <= 9; b += 2) {
      if ((a + b) % 3 == 0) {
        for (int test_d = 1; test_d <= MAX_DIGITS; test_d += 1) {
          // assert a * pow(10, test_d, 3) + b % 3 == 0
          is_prime[test_d][a][b] = 3;
        }
      }
    }
  }

  // Filter simple divisible by 7 mods.
  for (int test_d = 1; test_d <= MAX_DIGITS; test_d += 1) {
    is_prime[test_d][7][7] = 7;
  }

  for (int test_d = 3; test_d <= MAX_DIGITS; test_d += 1) {
    // See logic on Fermat primes:
    //   a^b + 1 can only be prime if b has no odd divisors
    //    => b is a power of two.
    bool is_power_two = (test_d & (test_d - 1)) == 0;
    if (!is_power_two) {
      is_prime[test_d][1][1] = -2; // Not prime but factor is unknown.
    }
  }
}

void FilterStats() {
  int total = 0;
  int total_to_test = 0;
  int filtered = 0;
  int filtered_trivial = 0;

  for (int d = START_DIGIT; d <= MAX_DIGITS; d++) {
    for (int a = 1; a <= 9; a++) {
      for (int b = 1; b <= 9; b++) {
        long status = is_prime[d][a][b];
        assert(status >= 0 || status == -2);

        total += 1;
        total_to_test += status == 0;
        filtered_trivial += (status >= 2 && status <= 5) || (a == 7 && b == 7);
        filtered += status != 0;
      }
    }
  }

  printf("%d total, %d trivially filtered, %d total filtered, %d to test (%.3f non-trivial remaining)\n",
         total, filtered_trivial, filtered, total_to_test, 1.0 * total_to_test / (total - filtered_trivial));
}

void VerifyFilter() {
  #pragma omp parallel for schedule( dynamic )
  for (long p = 2; p <= SIEVE_LIMIT; p++) {
    mpz_class p_mpz = p;
    if (mpz_probab_prime_p(p_mpz.get_mpz_t(), 25)) {
      mpz_class ten = 10;
      mpz_class t_mod;
      mpz_powm_ui(t_mod.get_mpz_t(), ten.get_mpz_t(), START_DIGIT - 1, p_mpz.get_mpz_t());

      long pow_ten_mod_p = t_mod.get_si();
      for (int d = START_DIGIT; d <= MAX_DIGITS; d++) {
        pow_ten_mod_p = (pow_ten_mod_p * 10) % p;
        if (pow_ten_mod_p == 0) {
          break;
        }

        for (long a = 1; a <= 9; a++) {
          for (long b = 1; b <= 9; b++) {
            // TODO: Deal with a * pow_ten_mod_p * b == p
            if (d <= 10) { continue; }

            long status = is_prime[d][a][b];
            if (status == 0) {
              if (((a * pow_ten_mod_p + b) % p) == 0) {
                cout << "ERROR: " << a << " * 10^" << d << " + " << b << " % " << p << " == 0" << endl;
              }
            }
          }
        }
      }
    }
  }
}

void SaveFilter() {
  char file_name[100];
  sprintf(file_name, "filter_%d_%d.filter", START_DIGIT, MAX_DIGITS);
  cout << "\tSaving to: " << file_name << endl;

  fstream fs;
  fs.open (file_name, std::fstream::out);

  // TODO: Record what prime divided filtered items.
  // TODO: same format as tester.cpp

  int count = 0;
  for (int d = START_DIGIT; d <= MAX_DIGITS; d++) {
    fs << d << ": ";
    for (long a = 1; a <= 9; a++) {
      for (long b = 1; b <= 9; b++) {
        long status = is_prime[d][a][b];
        if (status == 0 || d <= 10) {
          fs << "(" << a << "," << b << "), ";
          count += 1;
        }
      }
    }
    fs << endl;
  }

  fs.close();
  cout << endl;
  cout << "wc -w " << file_name << " - " << (MAX_DIGITS - START_DIGIT + 1)
       << " = " << count << " (number to test, includes extra small numbers)" << endl;
}
