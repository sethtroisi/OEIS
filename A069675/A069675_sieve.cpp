// Note the omp parallel for requires -fopenmp at compile time
// yields ~10x speedup.

#include <gmpxx.h>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

//#include <map>
//#include <google/dense_hash_map>
#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"

using namespace std;

#define START_DIGIT 1
#define MAX_DIGITS  200000

#define ONE_MILLION 1000000L
#define SIEVE_LIMIT 100 * ONE_MILLION

//                                   no d_step | d_step
// 40000, 1M   (147XXX to test):  Filter    ?s |
// 40000, 10M  (127069 to test):  Filter   13s |
// 40000, 20M  (121837 to test):  Filter   25s | 11.5s
// 40000, 100M (110499 to test):  Filter  111s | 52.5s
// 40000, 1B   (98337 to test):   Filter  982s |
// 40000, 2B   (95152 to test):   Filter 1899s |

// 200000, 10M    (630813 to test):   Filter    64s | 15s
// 200000, 100M   (551595 to test):   Filter      s | 44s
// 200000, 2B     (475227 to test):   Filter  9360s |

// ---- OLD ----
// NOTE the X to test are stale and ~50 higher because of a new test for a=b=1
// 3000, 10   (102501 to test): Filter  0s, User 8130, Parallel: 720
// 3000, 100  (31671 to test):  Filter  0s, User 4279, Parallel: 380
// 3000, 1000 (21654 to test):  Filter  0s, User 2998, Parallel: 263
// 3000, 100k (13303 to test):  Filter  0s, User 1914, Parallel: 167
// 3000, 1M   (11155 to test):  Filter  1s, User 1496, Parallel: 140
// 3000, 10M  (9550 to test):   Filter  7s, User 1481, Parallel: 135
// 3000, 20M  (9124 to test):   Filter 13s, User 1520, Parallel: 136
// 3000, 100M (8341 to test):   Filter 62s
// 3000, 1B   (7442 to test):   Filter 563s

// 5000, 10M  (15853 to test):  Filter  11s, User 8120, Parallel: 780
// 5000, 20M  (15175 to test):  Filter  30s, User 8087, Parallel: 732
// 5000, 30M  (14808 to test):  Filter  44s, User 8019, Parallel: 747
// 5000, 100M (13845 to test):  Filter 141s, User 8725, Parallel: 785

// 13000, 100M (36130 to test): Filter 231s, User 194031, Parallel: 16834

// 40000, 1M   (147509 to test):  Filter    9s,
// 40000, 10M  (126459 to test):  Filter   76s, | 24s
// 40000, 100M (110655 to test):  Filter  687s, | 203s
// 40000, 1B   (98483 to test):   Filter 5112s,
// 40000, 2B   (95248 to test):   Filter 9637s,

// 100000, 10B (12.... to test):  Filter 1440m

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

long is_prime[MAX_DIGITS+1][10][10] = {0};

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

void FilterSieve() {
  // Sieve out "small" prime factors and mark those numbers not to test.
  auto T0 = chrono::high_resolution_clock::now();

  vector<long> primes;
  primes.push_back(2);
  int prime_pi = 1;

  // Only odd indexes (divide by two to find value)
  vector<bool>test_p((SIEVE_LIMIT+1)/2, true);
  for (long p = 3; p*p <= SIEVE_LIMIT; p += 2) {
    if (test_p[p/2]) {
      for (long mi = p*p; mi <= SIEVE_LIMIT; mi += 2 * p) {
        test_p[mi/2] = false;
      }
    }
  }

  for (long p = 3; p <= SIEVE_LIMIT; p += 2) {
    if (test_p[p/2]) {
      primes.push_back(p);
      prime_pi += 1;
    }
  }

  auto T1 = chrono::high_resolution_clock::now();
  auto sieve_ms = chrono::duration_cast<chrono::milliseconds>(T1 - T0).count();
  cout << "PrimePi(" << ((SIEVE_LIMIT > ONE_MILLION) ? (SIEVE_LIMIT/ONE_MILLION) : SIEVE_LIMIT) <<
       ((SIEVE_LIMIT > ONE_MILLION) ? "M" : "") << ") = " << prime_pi <<
       " (" << sieve_ms << " ms)" << endl;

  // Instead of doing each power of ten, group multiple together
  // Do this by calculating a larger many divisible_mods at once
  // o(24 * d_step + d/dstep * hash_lookup(24 * d_step))
  // minimized = 24 - d / x^2 => x^2 = d / 24, x = sqrt(d/24)
  float ADJ_FACTOR = 0.6;
  int d_range = MAX_DIGITS - START_DIGIT + 1;
  int d_step = max(1, (int) (ADJ_FACTOR * sqrt(d_range / 24.0)));
  cout << "d_step of " << d_step << endl;

  mpz_class ten_d_step;
  mpz_ui_pow_ui(ten_d_step.get_mpz_t(), 10, d_step);

  #pragma omp parallel for schedule( dynamic )
  for (int pi = 0; pi < primes.size(); pi++) {
    long p = primes[pi];
    mpz_class p_mpz = p;

    if (pi * 20 % primes.size() < 20) { cout << "\tprime: " << p << endl; }

    if (p <= 5) {
      continue;
    }
    int count_divisible_mods = 0;
    //map<long, vector<tuple<int,int> > > divisible_mods;
    //
    //google::dense_hash_map<long, vector<tuple<int,int> > > divisible_mods;
    //divisible_mods.set_empty_key(-1);

    //absl::flat_hash_map<long, vector<tuple<int,int,int>>> divisible_mods;
    absl::flat_hash_set<long> divisible_mods;
    divisible_mods.reserve(24 * d_step);

    mpz_class ten_inverse_mpz;
    mpz_class ten_mpz = 10;
    mpz_invert(ten_inverse_mpz.get_mpz_t(), ten_mpz.get_mpz_t(), p_mpz.get_mpz_t());
    long ten_inverse = ten_inverse_mpz.get_si();

    for (long a = 1; a <= 9; a++) {
      if (a == p) {
        continue;
      }

      mpz_class inverse_a;
      mpz_class a_mpz = a;
      mpz_invert(inverse_a.get_mpz_t(), a_mpz.get_mpz_t(), p_mpz.get_mpz_t());

      // this is multiplied by 10 to cancel the first * ten_inverse (inside d loop).
      long inverse = (inverse_a.get_si() * 10) % p;
      for (int d = 0; d < d_step; d++) {
        inverse = (inverse * ten_inverse) % p;

        for (long b = 1; b <= 9; b += 2) {
          if ((b == 5) || ((a + b) % 3 == 0)) {
            continue;
          }

          // a * 10 ** d * p + b % p == 0
          // a * 10 ** h % p = - b
          long t = (inverse * -b) % p;

          //mpz_class t = (a_inverse * (p - b)) % p;
          //long t = (a_inverse.get_si() * -b) % p;
          if (t < 0) { t += p; };
          assert(t >= 0 && t < p);

          // (a * m_i) % p == 1
          // t = m_i * (p - b)
          //
          // a * t + b mod p =
          //   =  a * (m_i * (p - b)) + b
          //   =  a * m_i * p + a * m_i * -b + b
          //   =  1       * p + 1       * -b + b
          //   = -b + b
          //   = 0
          // =>
          //  (a * t + b) % p == 0
          //
          // if 10^d % p = t  =>  a * 10^d + b % p == 0
          //assert((a * t + b) % p == 0);
//          mpz_class ten_d_mod_mpz = 10;
//          mpz_powm_ui(ten_d_mod_mpz.get_mpz_t(), ten_d_mod_mpz.get_mpz_t(), d, p_mpz.get_mpz_t());
//          assert (((a * ten_d_mod_mpz * t + b) % p) == 0);

//          divisible_mods[t].push_back(make_tuple(d, a, b));
          divisible_mods.insert(t);
          count_divisible_mods += 1;
        }
      }
      //assert ((ten_inverse_mpz * inverse_ten_d * ten_d_step) % p == 1);
    }
    // if p is large, count_divisible_mods == 24 (9 * 4 * 2/3)
    assert (p < 100 || count_divisible_mods == 24 * d_step);

    int min_d = floor(log10(p));
    int start_d = max(min_d + 1, START_DIGIT);

    mpz_class ten_d_step_mod_mpz = ten_d_step % p;
    long ten_d_step_mod = ten_d_step_mod_mpz.get_si();
    long power_ten = 1;
    for (int d = 0; d <= MAX_DIGITS; d += d_step) {
      // THIS IS THE HOT BLOCK
      // Most of the computation happens here: o(primes * MAX_DIGITS/d_step) = O(billions)

      // Once had a optimization
      if (d > 0) {
        power_ten = (ten_d_step_mod * power_ten) % p;
        assert (power_ten > 0 && power_ten < p);
      }

      //if (d < START_DIGIT) { continue; }

      // This lookup takes 50-80% of all time.
      auto lookup = divisible_mods.find(power_ten);
      if (lookup != divisible_mods.end()) {
/*
        for (auto it = lookup->second.begin(); it != lookup->second.end(); it++)
          int d_inc = get<0>(*it);
          int a = get<1>(*it);
          int b = get<2>(*it);
*/
        // hits are relatively rare.
        for (long a = 1; a <= 9; a++) {
          long mod = a * power_ten % p;
          for (long d_inc = 0; d_inc < d_step; d_inc++) {
            if (d_inc != 0) {
              mod = (mod * 10) % p;
            }
            for (long b = 1; b <= 9; b += 2) {
              if (b == 5 || (a + b) % 3 == 0) { continue; }
              if ((mod + b) % p == 0) {
                int test_d = d + d_inc;
                if (test_d < START_DIGIT) { continue; }
                if (test_d > MAX_DIGITS) { continue; }

                if (is_prime[test_d][a][b] == 0) {
                  AssertDivisible(a, test_d, b, p);
                  is_prime[test_d][a][b] = p;
                }
              }
            }
          }
        }
      }
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
       << " = " << count << " (number to test)" << endl;
}


int main(void) {
  FilterSimple();

  auto T0 = chrono::high_resolution_clock::now();
  FilterSieve();
  FilterStats();
  // Takes 10-100x as FilterSieve()
  // VerifyFilter();

  SaveFilter();

  auto T1 = chrono::high_resolution_clock::now();
  auto filter_ms = chrono::duration_cast<chrono::milliseconds>(T1 - T0).count();
  cout << "Filter took " << filter_ms / 1000.0 << " seconds" << endl;
}
