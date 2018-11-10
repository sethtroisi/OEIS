// Note the omp parallel for requires -fopenmp at compile time
// yields ~4x speedup.

#include <gmpxx.h>
#include <chrono>
#include <cmath>
#include <iostream>

using namespace std;


int main() {
/*
167 3 * 10^67 + 1
218 7 * 10^263 + 9
253 7 * 10^783 + 9
277 6 * 10^1334 + 1
297 4 * 10^2245 + 3
310 9 * 10^2914 + 7
*/

// /*
  // Big primes
  vector<int> a = {3, 7, 7, 6, 4, 9};
  vector<int> d = {67, 263, 783, 1334, 2245, 2914};
  vector<int> b = {1, 9, 9, 1, 3, 7};
// */

/*
  vector <int> a = {2, 4, 6, 8};
  vector<int> d = {67, 263, 783, 1334, 2245, 2914};
  vector<int> b = {1, 0, 5, 1};
// */

  for (int reps = 0; reps < 30; reps++) {
    auto T0 = chrono::high_resolution_clock::now();

    int count = 0;
    for (int i = 0; i < a.size(); i++) {
      mpz_class t;
      mpz_ui_pow_ui(t.get_mpz_t(), 10, d[i]);
      t = a[i] * t + b[i];

      if (mpz_millerrabin(t.get_mpz_t(), reps)) {
        count += 1;
      }
    }

    auto T1 = chrono::high_resolution_clock::now();
    auto filter_ms = chrono::duration_cast<chrono::milliseconds>(T1 - T0).count(); 

    cout << "REPS: " << reps << " took " << filter_ms / 1000.0 << " seconds millerrabin, count: " << count << endl;


    T0 = chrono::high_resolution_clock::now();
    count = 0;
    for (int i = 0; i < a.size(); i++) {
      mpz_class t;
      mpz_ui_pow_ui(t.get_mpz_t(), 10, d[i]);
      t = a[i] * t + b[i];

      if (mpz_probab_prime_p(t.get_mpz_t(), reps)) {
        count += 1;
      }
    }

    T1 = chrono::high_resolution_clock::now();
    filter_ms = chrono::duration_cast<chrono::milliseconds>(T1 - T0).count(); 
    cout << "REPS: " << reps << " took " << filter_ms / 1000.0 << " seconds probab_prime_p, count: " << count << endl;
    cout << endl;
  }
}
