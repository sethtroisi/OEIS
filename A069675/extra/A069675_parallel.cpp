// Note the omp parallel for requires -fopenmp at compile time
// yields ~10x speedup.

#include <gmpxx.h>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

//#define START_DIGIT 1
//#define MAX_DIGITS 1000
#define START_DIGIT 35900
#define MAX_DIGITS  37500

#define ONE_MILLION 1000000
//#define SIEVE_LIMIT 1 * ONE_MILLION
#define SIEVE_LIMIT 1000 * ONE_MILLION

#define REPS 25

// N = 3000
// 0   => 2'30 = 150  (still does ~1 rep)
// 25  => 2'49 = 169
// 100 => 3'51 = 231

// Reps = 25
// 3000, 10   (102501 to test): Filter  0s, User 8130, Parallel: 720
// 3000, 100  (31671 to test):  Filter  0s, User 4279, Parallel: 380
// 3000, 1000 (21654 to test):  Filter  0s, User 2998, Parallel: 263
// 3000, 100k (13303 to test):  Filter  0s, User 1914, Parallel: 167
// 3000, 1M   (11155 to test):  Filter  1s, User 1496, Parallel: 140
// 3000, 10M  (9550 to test):   Filter  7s, User 1481, Parallel: 135
// 3000, 20M  (9124 to test):   Filter 13s, User 1520, Parallel: 136
//
// 5000, 10M  (15853 to test):  Filter  11s, User 8120, Parallel: 780
// 5000, 20M  (15175 to test):  Filter  30s, User 8087, Parallel: 732
// 5000, 30M  (14808 to test):  Filter  44s, User 8019, Parallel: 747
// 5000, 100M (13845 to test):  Filter 141s, User 8725, Parallel: 785


// 25, 13000, 100M (36130 to test): Filter 231s, User 194031, Parallel: 16834


// 25 24800-30000, 1B (12903 to test): Fliter 901s, User 1,247,280s Parallel: 108540


bool AssertDivisible(int a, int d, int b, long p) {
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


float GetCost(int a, int d, int b) {
  // ignores + b;
  float logN = log(a) + log(10) * d;

  // MillerRabin is O(log(n)^2)
  // In practice 2.3 makes eta more stable
  return pow(logN, 2.3);
}

long is_prime[MAX_DIGITS+1][10][10] = {0};
int approx_count = 4;

float TestD(int d) {
  mpz_class leading_one;
  mpz_ui_pow_ui(leading_one.get_mpz_t(), 10, d);

  float cost_finished = 0;

  for (int a = 1; a <= 9; a++) {
    mpz_class left = a * leading_one;

    for (int b = 1; b <= 9; b += 2) {
      if (is_prime[d][a][b] == 0) {
        assert(b != 5);
        assert((a + b) % 3 != 0);

        is_prime[d][a][b] = 1;
        cost_finished += GetCost(a, d, b);

        mpz_class t = left + b;

        if (mpz_millerrabin(t.get_mpz_t(), REPS)) {
        //if (mpz_probab_prime_p(t.get_mpz_t(), REPS)) {
          approx_count += 1;
          is_prime[d][a][b] = -1;

          // Small numbers are boring and clog the screen.
          if (d >= 3000) {
            cout << approx_count << " " << a << " * 10^" << d << " + " << b << endl;
          }
        }
      }
    }
  }
  return cost_finished;
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
    for (int b = 1; b <= 9; b += 1) {
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
  cout << "PrimePi(" << ((SIEVE_LIMIT > 2*ONE_MILLION) ? SIEVE_LIMIT/ONE_MILLION : SIEVE_LIMIT) <<
       ((SIEVE_LIMIT > ONE_MILLION) ? "M" : "") << ") = " << prime_pi <<
       " (" << sieve_ms << " ms)" << endl;

  #pragma omp parallel for schedule( dynamic )
  for (int pi = 0; pi < primes.size(); pi++) {
    long p = primes[pi];
    mpz_class p_mpz = p;

    if (p <= 5) {
      continue;
    }

    int count_divisible_mods = 0;
    map<long, vector<pair<int,int> > > divisible_mods;
    for (long a = 1; a <= 9; a++) {
      if (a == p) {
        continue;
      }

      mpz_class modular_inverse;
      mpz_class a_mpz = a;

      mpz_invert(modular_inverse.get_mpz_t(), a_mpz.get_mpz_t(), p_mpz.get_mpz_t());
      for (long b = 1; b <= 9; b += 2) {
        if ((b == 5) || ((a + b) % 3 == 0)) {
          continue;
        }

        mpz_class t = (modular_inverse * (p - b)) % p;
        if (t < 0) {
          t += p;
        }

        // a * t + b % p == 0
        // if 10^d % p = t  =>  a * 10^d + b % p == 0
        assert((a * t + b) % p == 0);
        divisible_mods[t.get_si()].push_back(make_pair(a, b));
        count_divisible_mods += 1;
      }
    }
    // if p is large, count_divisible_mods == 24 (9 * 4 * 2/3)


    // Technically we only need to go up to order but it doesn't save time.
    // Most primes have order > MAX_DIGITS

    int min_d = floor(log10(p));
    int start_d = max(min_d + 1, START_DIGIT);

    mpz_class t = 10;
    mpz_powm_ui(t.get_mpz_t(), t.get_mpz_t(), start_d - 1, p_mpz.get_mpz_t());

    long power_ten = t.get_si();
    for (int d = start_d; d <= MAX_DIGITS; d++) {
      power_ten = (10 * power_ten) % p;

      auto lookup = divisible_mods.find(power_ten);
      if (lookup != divisible_mods.end()) {
        for (auto it = lookup->second.begin(); it != lookup->second.end(); it++) {
          int a = it->first;
          int b = it->second;

          if (is_prime[d][a][b] == 0) {
            AssertDivisible(a, d, b, p);
            is_prime[d][a][b] = p;
          }
        }
      }
    }
  }
}

float FilterStats() {
  int total = 0;
  int total_to_test = 0;
  int filtered = 0;
  int filtered_trivial = 0;

  for (int d = START_DIGIT; d <= MAX_DIGITS; d++) {
    for (int a = 1; a <= 9; a++) {
      for (int b = 1; b <= 9; b++) {
        int status = is_prime[d][a][b];
        assert(status >= 0);

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

float PredictCost() {
  float total_cost = 0;

  for (int d = START_DIGIT; d <= MAX_DIGITS; d++) {
    for (int a = 1; a <= 9; a++) {
      for (int b = 1; b <= 9; b++) {
        int status = is_prime[d][a][b];
        assert(status >= 0);
        if (status == 0) {
          total_cost += GetCost(a, d, b);
        }
      }
    }
  }

  return total_cost;
}

void VerifyFilter() {
  #pragma omp parallel for schedule( dynamic )
  for (long p = 2; p <= SIEVE_LIMIT; p++) {
    mpz_class mpz_p = p;
    if (mpz_probab_prime_p(mpz_p.get_mpz_t(), REPS)) {
      mpz_class ten = 10;
      mpz_class t_mod;
      mpz_powm_ui(t_mod.get_mpz_t(), ten.get_mpz_t(), START_DIGIT - 1, mpz_p.get_mpz_t());

      long pow_ten_mod_p = t_mod.get_si();
      for (int d = START_DIGIT; d <= MAX_DIGITS; d++) {
        pow_ten_mod_p = (pow_ten_mod_p * 10) % p;
        if (pow_ten_mod_p == 0) {
          break;
        }

        for (long a = 1; a <= 9; a++) {
          for (long b = 1; b <= 9; b++) {
            // TODO deal with a * pow_ten_mod_p * b == p

            int status = is_prime[d][a][b];
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

int main(void) {
  cout << "Deprecated use A06975_ab_filter.cpp and A06975_tester.cpp" << endl;
  /*
  vector<string> seq;

  FilterSimple();

  float predicted_cost = 1.0 / 0.0; // positiy infinity

  auto T0 = chrono::high_resolution_clock::now();
  if (MAX_DIGITS > 400) {
    FilterSieve();
    FilterStats();
    // Takes 10-100x as FilterSieve()
    //VerifyFilter();

    predicted_cost = PredictCost();
    auto T1 = chrono::high_resolution_clock::now();
    auto filter_ms = chrono::duration_cast<chrono::milliseconds>(T1 - T0).count();
    cout << "\tFilter took " << filter_ms / 1000.0 << " seconds" << endl;
    if (START_DIGIT > 1) {
      printf("\tO(%d) - O(%d) ~= %.2f\n", MAX_DIGITS, START_DIGIT, predicted_cost);
    } else {
      printf("\tO(%d) ~= %.2f\n", MAX_DIGITS, predicted_cost);
    }
    cout << endl;
  }
  auto T2 = chrono::high_resolution_clock::now();

  float cost_done = 0;

  // Deal with fast (and numerous initial values)
  for (int d = START_DIGIT; d <= 200; d++) {
    cost_done += TestD(d);
  }

  // Today make predicted_cost return approx seconds and print once every 10
  // minutes.
  float status_prints = 10;
  if (MAX_DIGITS > 1000) { status_prints = 50; }
  if (MAX_DIGITS > 5000) { status_prints = 100; }
  if (MAX_DIGITS > 10000) { status_prints = 500; }
  float cost_done_print = predicted_cost / status_prints;

  #pragma omp parallel for schedule( dynamic )
  for (int d = max(201, START_DIGIT); d <= MAX_DIGITS; d++) {
    cost_done += TestD(d);

    if (cost_done >= cost_done_print) {
      auto T3 = chrono::high_resolution_clock::now();
      auto elapsed_ms = chrono::duration_cast<chrono::milliseconds>(T3 - T2).count();
      float eta_ms = (predicted_cost / cost_done) * elapsed_ms;

      printf("Finished d: %d,  %.1f%% (%ld seconds, estimate: %.0f seconds)\n",
             d, 100 * cost_done / predicted_cost, elapsed_ms / 1000, eta_ms / 1000);
      cost_done_print += predicted_cost / status_prints;
    }
  }

  auto T4 = chrono::high_resolution_clock::now();
  auto primality_ms = chrono::duration_cast<chrono::milliseconds>(T4 - T2).count();
  cout << "\tCheck primality took " << primality_ms / 1000.0 << " seconds" << endl;
  cout << endl << endl;


  int count = 0;
  if (START_DIGIT == 1) {
    for (mpz_class t = 1; t < 10; t++) {
      if (mpz_probab_prime_p(t.get_mpz_t(), REPS)) {
        count += 1;
        //seq.push_back(t.get_str());

        //cout << count << " " << t << endl;
      }
    }
  }
  for (int d = START_DIGIT; d <= MAX_DIGITS; d++) {
    for (int a = 1; a <= 9; a++) {
      for (int b = 1; b <= 9; b += 2) {
        if (is_prime[d][a][b] == 0) {
          cout << "PROBLEM: " << d << " " << a << " " << b << endl;
        }

        if (is_prime[d][a][b] == -1) {
          count += 1;
          //mpz_class t = a * leading_one + b;
          //seq.push_back(t.get_str());

          if (count % 10 == 0 || d >= 3000) {
            cout << count << " " << a << " * 10^" << d << " + " << b << endl;
          }
        }
      }
    }
  }
  //bfile.WriteListToFile("069675", seq)
  */
}
