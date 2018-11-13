// Note the omp parallel for requires -fopenmp at compile time
// yields ~10x speedup.

#include <gmpxx.h>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "A069675_gpu_shared.h"
#include "A069675_gpu_cuda_sieve.h"

using namespace std;

// NOTE the X to test are stale and ~50 higher because of a new test for a=b=1

// 3000, 100  (31671 to test):  Filter  0s, User 4279, Parallel: 380
// 3000, 1M   (11155 to test):  Filter  1s, User 1496, Parallel: 140
// 3000, 10M  (9550 to test):   Filter  7s, User 1481, Parallel: 135
// 3000, 100M (8341 to test):   Filter 62s
// 3000, 1B   (7442 to test):   Filter 563s

// 243000 total, 174000 trivially filtered, 177496 total filtered, 65504 to test (0.949 non-trivial remaining)


// 5000, 10M  (15853 to test):  Filter  11s, User 8120, Parallel: 780
// 5000, 100M (13845 to test):  Filter 141s, User 8725, Parallel: 785

// 40000, 1M   (147509 to test):  Filter    9s,
// 40000, 10M  (126459 to test):  Filter   76s, | 24s
// 40000, 100M (110655 to test):  Filter  687s, | 203s
// 40000, 1B   (98483 to test):   Filter 5112s,
// 40000, 2B   (95248 to test):   Filter 9637s,

// 100000, 10B (12.... to test):  Filter 1440m

// 405000 total, 290000 trivially filtered, 295829 total filtered, 109171 to test (0.949 non-trivial remaining)


/*
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
*/

long is_prime[MAX_DIGITS_P1][10][10] = {0};

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
  {
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
    assert(prime_pi = primes.size());
  }

  auto T1 = chrono::high_resolution_clock::now();
  auto sieve_ms = chrono::duration_cast<chrono::milliseconds>(T1 - T0).count();
  cout << "PrimePi(" << ((SIEVE_LIMIT > ONE_MILLION) ? (SIEVE_LIMIT/ONE_MILLION) : SIEVE_LIMIT)
       << ((SIEVE_LIMIT > ONE_MILLION) ? "M" : "") << ") = "
       << prime_pi << " (" << sieve_ms << " ms)" << endl;

  // NOTE: This will be the first thing that needs to be segmentede
  // TODO make sure this works with big prime_pi
  auto div_mods = new long[prime_pi][24];
  if (div_mods == nullptr) { exit(1); }

  #pragma omp parallel for schedule( dynamic )
  for (int pi = 0; pi < prime_pi; pi++) {
    long p = primes[pi];

    if (p <= 5) {
      continue;
    }

    int count_divisible_mods = 0;
    for (long a = 1; a <= 9; a++) {
      if (a == p) {
        continue;
      }

      mpz_class modular_inverse;
      mpz_class a_mpz = a;
      mpz_class p_mpz = p;
      mpz_invert(modular_inverse.get_mpz_t(), a_mpz.get_mpz_t(), p_mpz.get_mpz_t());

      for (long b = 1; b <= 9; b += 2) {
        if ((b == 5) || ((a + b) % 3 == 0)) {
          continue;
        }

        // NOTE: Math in A069675_sieve.cpp
        long t = (modular_inverse.get_si() * -b) % p;
        if (t < 0) { t += p; };
        assert(t >= 0 && t < p);

        div_mods[pi][count_divisible_mods] = t;
        count_divisible_mods += 1;

        // These can be 'recovered' by just testing a * t + b for all a,b
        // divisible_mods[t].push_back(make_pair(a, b));
      }
    }
    sort(div_mods[pi], div_mods[pi] + count_divisible_mods);

    // if p is large, count_divisible_mods == 24 (9 * 4 * 2/3)
    assert (p < 10 ? (count_divisible_mods <= 24) : (count_divisible_mods == 24));
  }

  auto T2 = chrono::high_resolution_clock::now();
  auto filter_pre_ms = chrono::duration_cast<chrono::milliseconds>(T2 - T1).count();
  cout << "Filter pre-work " << filter_pre_ms << " ms" << endl;

  // Do small primes on host first
  int small_pi_limit = min(prime_pi, PRIME_PI_1M);
  #pragma omp parallel for schedule( dynamic )
  for (int pi = 0; pi < small_pi_limit; pi++) {
    if (primes[pi] <= 5) {
      continue;
    }

    test_p((long*)is_prime, primes[pi], div_mods[pi]);
  }

  // TODO turn this into a #define.
  auto T3 = chrono::high_resolution_clock::now();
  auto small_test_p_ms = chrono::duration_cast<chrono::milliseconds>(T3 - T2).count();
  cout << "Small primes " << small_test_p_ms << " ms" << endl;

  if (prime_pi > small_pi_limit) {
    auto results = new bool[prime_pi];
    memset(results, 0, prime_pi);

    FilterSieveKernelHost(
        is_prime,
        div_mods,
        primes.data(),
        small_pi_limit,
        prime_pi,
        results);

    auto T4 = chrono::high_resolution_clock::now();
    auto gpu_ms = chrono::duration_cast<chrono::milliseconds>(T4 - T3).count();
    cout << "GPU " << gpu_ms << " ms" << endl;

    long a = 0, b = 0;
    for (int pi = small_pi_limit + 1; pi < prime_pi; pi++) {
      a += 1;
      if (results[pi]) {
        b += 1;
        test_p((long*)is_prime, primes[pi], div_mods[pi]);
      }
    }
    if (a > 0) {
      printf("GPU filtered %ld to %ld (%.3f)\n", a, b, 1.0 * a / b);
    }

    delete[] results;
  }

  delete[] div_mods;
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
        if (status == 0) {
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

  //SaveFilter();

  auto T1 = chrono::high_resolution_clock::now();
  auto filter_ms = chrono::duration_cast<chrono::milliseconds>(T1 - T0).count();
  cout << "Filter took " << filter_ms / 1000.0 << " seconds" << endl;
}
