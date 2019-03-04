// Note the omp parallel for requires -fopenmp at compile time
// yields ~10x speedup.

#include <gmpxx.h>

#include <algorithm>
#include <atomic>
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

#include "A069675_sieve.h"
#include "A069675_extra.h"

// defensive programming makes everything 10-20% slower an acceptable
// reduction. This can be modified by toggling the cmake file.

// If this is not true than mod * mod may overflow int64
// Makes everythings ~2.5 slower (int64 replaced with mpz_class)
#define IS_SMALL_PRIMES (SIEVE_LIMIT < (2001 * ONE_MILLION))
//#define IS_SMALL_PRIMES false

// 0.55 with no large prime support
// 0.65 with large prime support
#define ADJ_FACTOR 0.75

#define SEGMENTS_SIZE (500L * ONE_MILLION)

#define CKPT_PER_SEGMENT 5

using namespace std;

//                                   no d_step | d_step
// 40000, 1M   (147XXX to test):  Filter    ?s |
// 40000, 10M  (127069 to test):  Filter   13s |
// 40000, 20M  (121837 to test):  Filter   25s | 11.5s
// 40000, 100M (110499 to test):  Filter  111s | 52.5s
// 40000, 1B   (98337 to test):   Filter  982s |
// 40000, 2B   (95152 to test):   Filter 1899s |

// 200000, 10M    (630813 to test):   Filter    64s | 5.3s
// 200000, 50M    (573124 to test):   Filter        | 22s  (ADJ 0.55 small)    48s (0.75)
// 200000, 100M   (551595 to test):   Filter        | 39s  (ADJ 0.55 small)    74s (0.75)
// 200000, 1B     (490438 to test):   Filter        | 334s (ADJ 0.55 small)   746s (0.75)
// 200000, 2B     (474425 to test):   Filter  9360s | 646s

// ---- OLD ----
// NOTE the X to test are stale and ~50 higher because of a new test for a=b=1
// 5000, 10M  (15853 to test):  Filter  11s, User 8120, Parallel: 780
// 5000, 100M (13845 to test):  Filter 141s, User 8725, Parallel: 785

// 40000, 1M   (147509 to test):  Filter    9s,
// 40000, 10M  (126459 to test):  Filter   76s,
// 40000, 100M (110655 to test):  Filter  687s,
// 40000, 1B   (98483 to test):   Filter 5112s,
// 40000, 2B   (95248 to test):   Filter 9637s,

// 100000, 10B (12.... to test):  Filter 1440m

long is_prime[MAX_DIGITS+1][10][10] = {0};
atomic<int> sieve_filtered(0);

void filterP(long p, const long d_step, const mpz_class& ten_d_step) {
  mpz_class p_mpz = p;

  if (p <= 5) {
    // Handled by filter simple.
    return;
  }

  // Maps t to (a,b) pair
  // Tested array, map, google::dense_hash_map, absl::flat_hash_map
  absl::flat_hash_map<long, vector<tuple<int,int>>> divisible_mods_d1;
  divisible_mods_d1.reserve(24);

  // keys from divisible_mods for d = 0 to d_step (multiplied by inverse_ten)
  absl::flat_hash_set<long> divisible_mods;
  divisible_mods.reserve(24 * d_step);
  int count_divisible_mods = 0;

  mpz_class ten_inverse_mpz;
  mpz_class ten_mpz = 10;
  mpz_invert(ten_inverse_mpz.get_mpz_t(), ten_mpz.get_mpz_t(), p_mpz.get_mpz_t());
  long ten_inverse = ten_inverse_mpz.get_si();

  for (long a = 1; a <= 9; a++) {
    if (a == p) {
      continue;
    }

    mpz_class inverse_a_mpz;
    mpz_class a_mpz = a;
    mpz_invert(inverse_a_mpz.get_mpz_t(), a_mpz.get_mpz_t(), p_mpz.get_mpz_t());
    long inverse_a = inverse_a_mpz.get_si();

    // save t such that (d,a,b) indicates factor of p
    for (long b = 1; b <= 9; b += 2) {
      if ((b == 5) || ((a + b) % 3 == 0)) {
        continue;
      }

      long t = (inverse_a * -b) % p;
      if (t < 0) { t += p; };
      assert(t >= 0 && t < p);

      // (a * i_a) % p == 1
      // =>
      //  (a * t + b) % p == 0
      // if 10^d % p = t
      // =>
      //  a * 10^d + b % p == 0

      assert((a * t + b) % p == 0);

      divisible_mods_d1[t].push_back(make_tuple(a, b));
    }

    // save t such that (d,a,b) indicates factor of p
    // this is multiplied by 10 to cancel the first * ten_inverse (inside d loop).
    mpz_class inverse_mpz = (inverse_a_mpz * 10) % p;
    long inverse = (inverse_a * 10) % p;
    long ten_d_mod_p = 1;
    for (int d = 0; d < d_step; d++) {
      #if IS_SMALL_PRIMES
        inverse = (inverse * ten_inverse) % p;
      #else
        inverse_mpz = (inverse_mpz * ten_inverse_mpz) % p;
        inverse = inverse_mpz.get_si();
      #endif
      assert (inverse > 0 && inverse < p);

      for (long b = 1; b <= 9; b += 2) {
        if ((b == 5) || ((a + b) % 3 == 0)) {
          continue;
        }

        long t = (inverse * -b) % p;
        if (t < 0) { t += p; };

        // This is probably the most expensive assert
        // it gets enabled / disabled by hand
        // assert((((a * ten_d_mod_p) % p) * t + b) % p == 0);
        divisible_mods.insert(t);
        count_divisible_mods += 1;
      }

      ten_d_mod_p = (ten_d_mod_p * 10) % p;
    }
  }
  // if p is large, count_divisible_mods == 24 == 9 * 4 * 2/3
  assert (p < 100 || count_divisible_mods == 24 * d_step);

  // SETUP COMPLETE

  int min_d = floor(log10(p));
  int start_d = max(min_d + 1, START_DIGIT);

  // If d's that matched divisible_mods_d1 were cached,
  // loop could be exited early when cycle (e.g order) was found
  // Would only have 50% work (see Mathematica Notebook)

  mpz_class ten_d_step_mod_mpz = ten_d_step % p;
  long ten_d_step_mod = ten_d_step_mod_mpz.get_si();
  mpz_class power_ten_mod_p_mpz = 1;
  long power_ten_mod_p = 1;
  for (int d = 0; d <= MAX_DIGITS; d += d_step) {
    if (d != 0) {
      #if IS_SMALL_PRIMES
        power_ten_mod_p = (ten_d_step_mod * power_ten_mod_p) % p;
      #else
        // MUCH (4x?) slower than in64 multiple but supports p > 2B.
        // I suggest tweaking up ADJ_FACTOR
        power_ten_mod_p_mpz = (ten_d_step_mod * power_ten_mod_p_mpz) % p;
        power_ten_mod_p = power_ten_mod_p_mpz.get_si();
      #endif
      assert (power_ten_mod_p > 0 && power_ten_mod_p < p);
    }

    //if (d < START_DIGIT) { continue; }

    // This lookup takes 50-80% of all time.
    // And hits are relatively rare (after small p).
    auto lookup = divisible_mods.find(power_ten_mod_p);
    if (lookup != divisible_mods.end()) {
      // Something in d = 0 to d_step, a = 1 to 9, b odd will divide by p
      // Tried stored all (d,a,b) with divisible_mods
      // Tried scanning all (d,a,b) range
      // Now scanning d range with (a,b) lookup, seems fastest

      long temp = power_ten_mod_p;
      for (long d_inc = 0; d_inc < d_step; d_inc++) {
        if (d_inc != 0) {
          temp = (temp * 10) % p;
        }
        int test_d = d + d_inc;
        if (test_d < START_DIGIT) { continue; }
        if (test_d > MAX_DIGITS) { continue; }

        auto lookup_ab = divisible_mods_d1.find(temp);
        if (lookup_ab != divisible_mods_d1.end()) {
          for (const auto &it : lookup_ab->second) {
            int a = get<0>(it);
            int b = get<1>(it);
            if (is_prime[test_d][a][b] == 0) {
              AssertDivisible(a, test_d, b, p);
              is_prime[test_d][a][b] = p;
              sieve_filtered++;
            }
          }
        }
      }
    }
  }
}


vector<long> SegmentedSieveOfErat(long start, long stop, const vector<long>& small_primes) {
  // inclusive of start and stop
  auto T0 = chrono::high_resolution_clock::now();

  vector<long> primes;
  if (start <= 2 && stop >= 2) {
    primes.push_back(2);
  }

  // first odd after start
  long odd_start = start + (start % 2 == 0);
  assert(odd_start % 2 == 1);

  // last odd before end;
  long odd_stop = stop - (stop % 2 == 0);
  assert(odd_stop % 2 == 1);

  // Only odd values.
  // odd n is stored at (n - odd_start) / 2
  long indexes = (odd_stop - odd_start) / 2 + 1;
  vector<bool>test_p(indexes, true);

  if (odd_start == 1) {
    // 0 maps to 1 which is not prime
    test_p[0] = false;
  }

  // Skip 2
  for (long p : small_primes) {
    if (p == 2) { continue; }

    long first = p * p;
    if (first > stop) { break; }

    if (first < start) {
      // first odd multiple of p >= start
      long mult = (start - 1) / p + 1;
      mult += 1 - (mult % 2);

      first = p * mult;

      assert(first % 2 == 1);
      assert(first >= start);
      assert(first - 2*p < start);
    }

    // Multiple index
    long mi = (first - odd_start) / 2;
    for (; mi < indexes; mi += p) {
      test_p[mi] = false;
    }
  }

  // first odd after start
  for (long i = 0; i < indexes; i++) {
    if (test_p[i]) {
      long p = odd_start + 2 * i;
      primes.push_back(p);
      assert(start <= p && p <= stop);
    }
  }

  auto T1 = chrono::high_resolution_clock::now();
  auto sieve_ms = chrono::duration_cast<chrono::milliseconds>(T1 - T0).count();

  if ((start % ONE_BILLION == 0) && (stop % ONE_BILLION == 0)) {
    cout << "PrimePi("
         << start / ONE_BILLION << "B, "
         << stop  / ONE_BILLION << "B) = ";
  } else {
    cout << "PrimePi(" << start << ", " << stop << ") = ";
  }
  cout << primes.size()
       << "  (" << sieve_ms << " ms)" << endl;

  return primes;
}

vector<long> SmallSieveOfErat(long stop) {
  vector<long> primes;
  primes.push_back(2);

  // Only odd values.
  // odd n is stored at n / 2
  long indexes = (stop + 1) / 2;
  vector<bool>test_p(indexes, true);

  for (long p = 3; p*p <= stop; p += 2) {
    if (test_p[p/2]) {
      for (long mi = p*p; mi <= stop; mi += 2 * p) {
        test_p[mi/2] = false;
      }
    }
  }

  for (long p = 3; p <= stop; p += 2) {
    if (test_p[p/2]) {
      primes.push_back(p);
    }
  }
  return primes;
}

void FilterSieve() {
  auto T0 = chrono::high_resolution_clock::now();

  // Instead of doing each power of ten, group multiple together
  // Do this by calculating a larger many divisible_mods at once
  // o(24 * d_step + d/dstep * hash_lookup(24 * d_step))
  // minimized = 24 - d / x^2 => x^2 = d / 24, x = sqrt(d/24)
  int d_range = MAX_DIGITS - START_DIGIT + 1;
  int d_step = max(1, (int) (ADJ_FACTOR * sqrt(d_range / 24.0)));
  cout << endl;
  cout << "\tUsing d_step = " << d_step << endl << endl;

  mpz_class ten_d_step;
  mpz_ui_pow_ui(ten_d_step.get_mpz_t(), 10, d_step);

  // Sieve out "small" prime factors and mark those numbers not to test.
  auto sieve_primes = SmallSieveOfErat(ONE_MILLION);

  long start = 0;
  for (long stop = SEGMENTS_SIZE; start < SIEVE_LIMIT; stop += SEGMENTS_SIZE) {
    stop = min(stop, SIEVE_LIMIT);

    // Segmented sieve.
    auto primes = SegmentedSieveOfErat(start, stop, sieve_primes);

    #pragma omp parallel for schedule( dynamic )
    for (int pi = 0; pi < primes.size(); pi++) {
      long p = primes[pi];
      if (pi * CKPT_PER_SEGMENT % primes.size() < CKPT_PER_SEGMENT) {
        auto T1 = chrono::high_resolution_clock::now();
        auto sieve_s = chrono::duration_cast<chrono::seconds>(T1 - T0).count();
        cout << "\tprime(" << pi << ") = " << p << "  @" << sieve_s << endl;
      }

      filterP(p, d_step, ten_d_step);
    }

    cout << "\tfiltered " << sieve_filtered << " from primes <= " << stop << endl;

    // double counts start but it's even so it doesn't add a prime.
    start = stop;

    // Partial Status saving. COuld be make more frequent or something
    SaveFilter();
  }
  auto T1 = chrono::high_resolution_clock::now();
  auto filter_ms = chrono::duration_cast<chrono::milliseconds>(T1 - T0).count();
  cout << "Filter(" << SIEVE_LIMIT << ") took " << filter_ms / 1000.0 << " seconds" << endl;
}

int main(void) {
  cout << endl;
  #if !IS_SMALL_PRIMES
    cout << "Large Prime support needed, ~2x performance loss :(" << endl;
  #endif
  cout << endl;


  sieve_filtered = FilterSimple();
  cout << "Filtered " << sieve_filtered << " trivially" << endl;

  auto T0 = chrono::high_resolution_clock::now();

  FilterSieve();

  FilterStats();
  VerifyFilter();

  SaveFilter();

  auto T1 = chrono::high_resolution_clock::now();
  auto total_ms = chrono::duration_cast<chrono::milliseconds>(T1 - T0).count();
  cout << "Program took " << total_ms / 1000.0 << " seconds" << endl;
}
