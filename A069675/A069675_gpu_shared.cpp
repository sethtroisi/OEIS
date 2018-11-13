// Note the omp parallel for requires -fopenmp at compile time
// yields ~10x speedup.

#include "A069675_gpu_shared.h"

#include <array>
#include <cassert>
#include <cstdio>
#include <vector>

using namespace std;

CUDAHOSTDEV bool test_p(
    long is_prime_ref[MAX_DIGITS][10][10],
    long p,
    long mods[24]) {
  // Start at d = 20, everything smaller is fast
  int start_d = max(20, START_DIGIT);

  long power_ten = 1;
  for (int d = 1; d <= MAX_DIGITS; d++) {
    // THIS IS THE HOT BLOCK
    // Most of the computation happens here: o(primes * MAX_DIGITS) = O(billions)

    power_ten = (10 * power_ten) % p;

    if (d < START_DIGIT) { continue; }

    // In C++ with hashtable this takes 50-80% of all time.
    // TODO(seth): Investigate hand coded binary search or skipping by two.
    for (int i = 0; i < 24 && mods[i] <= power_ten; i++) {
      if (mods[i] == power_ten) {
        for (long a = 1; a <= 9; a++) {
          for (long b = 1; b <= 9; b++) {
            if ((a * power_ten + b) % p == 0) {
              if (is_prime_ref[d][a][b] == 0) {
                is_prime_ref[d][a][b] = p;
                #ifdef __CUDA__ARCH__
                  // On CUDA exit early and let host do full sweep.
                  return true;
//                #else
//                  AssertDivisible
                #endif
              }
            }
          }
        }
      }
    }
  }
  return false;
}


    // Avoiding division by doing max 3 subtrations helped in C++
    /*
    power_ten = 10 * power_ten;
    // 0 <= power_ten <= 10 * p, try to subtract 8p, 4p, 2p, 1p
    long shift_p = p << 3;

    // This should be faster but is 20-40% worse, why?
    // while (power_ten >= p) {
    while (shift_p >= p) {
      if (power_ten >= shift_p) {
        power_ten -= shift_p;
      }
      shift_p >>= 1;
    }
    */
    //assert (power_ten > 0 && power_ten < p);

