// Note the omp parallel for requires -fopenmp at compile time
// yields ~10x speedup.

#include <gmpxx.h>

#include <atomic>
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

int FilterSimple() {
  long filtered = 0;

  // Filter divisible by 2 and 5 mods
  for (int a = 1; a <= 9; a += 1) {
    for (int test_d = 1; test_d <= MAX_DIGITS; test_d += 1) {
      is_prime[test_d][a][2] = 2;
      is_prime[test_d][a][4] = 2;
      is_prime[test_d][a][6] = 2;
      is_prime[test_d][a][8] = 2;
      is_prime[test_d][a][5] = 5;
      filtered += 5;
    }
  }

  // Filter divisible by 3 mods.
  for (int a = 1; a <= 9; a += 1) {
    for (int b = 1; b <= 9; b += 2) {
      if ((a + b) % 3 == 0) {
        for (int test_d = 1; test_d <= MAX_DIGITS; test_d += 1) {
          // assert a * pow(10, test_d, 3) + b % 3 == 0
          if (is_prime[test_d][a][b] == 0) {
            is_prime[test_d][a][b] = 3;
            filtered++;
          }
        }
      }
    }
  }

  {
    // Filter simple, divisible by 7.
    // 10 ** d % 7 = 1,3,2,6,4,5, (repeats)
    int ten_d_mod_seven[] = {1,3,2,6,4,5};
    for (int test_d = 1; test_d <= MAX_DIGITS; test_d += 1) {
      int d_mod = ten_d_mod_seven[test_d % 6];
      for (int a = 1; a <= 9; a += 1) {
        for (int b = 1; b <= 9; b += 2) {
          if ((a * d_mod + b) % 7 == 0 && is_prime[test_d][a][b] == 0) {
            is_prime[test_d][a][b] = 7;
            filtered++;
          }
        }
      }
    }
  }

  for (int test_d = 3; test_d <= MAX_DIGITS; test_d += 1) {
    // See logic on Fermat primes:
    //   a^b + 1 can only be prime if b has no odd divisors
    //    => b is a power of two.
    bool is_power_two = (test_d & (test_d - 1)) == 0;
    if (!is_power_two) {
      is_prime[test_d][1][1] = -2; // Not prime but factor is unknown.
      filtered++;
    }
  }
  return filtered;
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

  cout << endl;
  printf("%d total, %d trivially filtered, %d total filtered\n",
         total, filtered_trivial, filtered);
  printf("\t%d to test (%.3f non-trivial remaining)\n",
         total_to_test, 1.0 * total_to_test / (total - filtered_trivial));
}

void VerifyFilter() {
  atomic<long> verified(0);
  atomic<long> negative(0);

  // No longer tests that things weren't missed, only that all entries divide.
  #pragma omp parallel for schedule( dynamic )
  for (int d = START_DIGIT; d <= MAX_DIGITS; d++) {
    if (d <= 10) { continue; }

    long v = 0;
    long n = 0;

    mpz_class ten = 10;
    mpz_pow_ui(ten.get_mpz_t(), ten.get_mpz_t(), d);

    for (long a = 1; a <= 9; a++) {
      for (long b = 1; b <= 9; b++) {
        // TODO: Deal with a * pow_ten_mod_p * b == p

        long status = is_prime[d][a][b];
        long p = status;
        if (status > 0) {
          v += 1;
          mpz_class modulo = (a * ten + b) % p;
          if (modulo != 0) {
            cout << "ERROR: " << a << " * 10^" << d << " + " << b << " % " << p
                 << " == " << modulo << endl;
          }
        } else if (status < 0) {
          n += 1;
        }
      }
    }
    verified += v;
    negative += n;
  }

  cout << endl;
  printf("verified %ld entries, skipped %ld negative entries\n",
         verified.load(), negative.load());
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
  cout << "wc -w " << file_name << " - " << (MAX_DIGITS - START_DIGIT + 1)
       << " = " << count << " (number to test, includes extra small numbers)" << endl;
  cout << endl;
}
