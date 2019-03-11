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

using namespace std;

#define START_DIGIT 1
#define MAX_DIGITS  200000

#define REPS 25

// N = 3000
// Reps = 0   => 2'30 = 150  (still does ~1 rep)
// Reps = 25  => 2'49 = 169
// Reps = 100 => 3'51 = 231

// Reps = 25
// 3000, 10   (102501 to test): Filter  0s, User 8130, Parallel: 720
// 3000, 100  (31671 to test):  Filter  0s, User 4279, Parallel: 380
//           | 30978 |                  0s,                           230
// 3000, 1000 (21654 to test):  Filter  0s, User 2998, Parallel: 263
//           | 21595 |                  0s,                           155
// 3000, 100k (13303 to test):  Filter  0s, User 1914, Parallel: 167
// 3000, 1M   (11155 to test):  Filter  1s, User 1496, Parallel: 140
//           | 11145 |                  0s,                           85 (+40%)
// 3000, 10M  (9550 to test):   Filter  7s, User 1481, Parallel: 135
// 3000, 20M  (9124 to test):   Filter 13s, User 1520, Parallel: 136
//           | 9118  |                  4s,                           72

// 5000, 10M  (15853 to test):  Filter  11s, User 8120, Parallel: 780
//           | 15839 |                  3s,                           455
// 5000, 20M  (15175 to test):  Filter  30s, User 8087, Parallel: 732
//           | 15166 |                  7s,                           415
// 5000, 30M  (14808 to test):  Filter  44s, User 8019, Parallel: 747
//           | 14799 |                  9s,                           404
// 5000, 100M (13845 to test):  Filter 141s, User 8725, Parallel: 785
//           | 13837 |                 30s,                           383

// 25, 13000, 100M (36130 to test): Filter 231s, User 194031, Parallel: 16834
//           | 36103 |                 68s,                           10145

// 25, 40000, 2B


float GetCost(int a, int d, int b) {
  // ignores + b;
  float logN = log(a) + log(10) * d;

  // MillerRabin is O(log(n)^2)
  // In practice 2.3 makes eta more stable
  return pow(logN, 2.3);
}

bool new_status[MAX_DIGITS+1][10][10] = {};
long is_prime[MAX_DIGITS+1][10][10] = {};
int approx_count = 4;

float TestD(int d) {
  mpz_class leading_one;
  mpz_ui_pow_ui(leading_one.get_mpz_t(), 10, d);

  float cost_finished = 0;

  for (int a = 1; a <= 9; a++) {
    mpz_class left = a * leading_one;

    // TODO consider adding parell here instead of outer (so that each TestD
    // finisheds faster and potentially more in order)
    for (int b = 1; b <= 9; b += 2) {
      if (is_prime[d][a][b] == 0) {
        //assert(b != 5);
        //assert((a + b) % 3 != 0);

        mpz_class t = left + b;

        if (mpz_millerrabin(t.get_mpz_t(), REPS)) {
        //if (mpz_probab_prime_p(t.get_mpz_t(), REPS)) {
          approx_count += 1;
          is_prime[d][a][b] = -1; // Prime!

          // Small numbers are boring and clog the screen.
          if (d >= 1000) {
            cout << approx_count << " " << a << " * 10^" << d << " + " << b << endl;
          }
        } else {
          is_prime[d][a][b] = -2; // Composite but no idea of factor.
        }
        new_status[d][a][b] = true;
        cost_finished += GetCost(a, d, b);
      }
    }
  }
  return cost_finished;
}

void LoadFilter() {
  string file_name = "filter_" + to_string(START_DIGIT) + "_" + to_string(MAX_DIGITS) + ".filter";
  cout << "\t\tReading from: " << file_name << endl;
  FILE *fs = fopen(file_name.c_str(), "r");

  for (int d = START_DIGIT; d <= MAX_DIGITS; d++) {
    string line;
    int testD;
    assert(fscanf(fs, "%d: ", &testD) == 1);
    assert (testD == d);

    for (long a = 1; a <= 9; a++) {
      for (long b = 1; b <= 9; b++) {
        is_prime[d][a][b] = 1; // Don't test
      }
    }

    int a, b;
    while (2 == fscanf(fs, "(%d, %d), ", &a, &b)) {
      is_prime[d][a][b] = 0; // Test these
      //cout << "\t" << a << ", " << b << endl;
    }
  }
  fclose(fs);
}

void LoadPartialResults() {
  string file_name = "filter_" + to_string(START_DIGIT) + "_" + to_string(MAX_DIGITS) + ".partial";
  FILE *fs = fopen(file_name.c_str(), "r");
  if (!fs) {
    cout << "\t\tno partial file @ " << file_name << endl;
    return;
  }
  cout << "\t\tReading from: " << file_name << endl;

  int testD, testA, testB;
  long testResult;
  while (4 == fscanf(fs, "%d, %d, %d: %ld", &testD, &testA, &testB, &testResult)) {
    assert(testResult <= -1);
    is_prime[testD][testA][testB] = testResult;
  }
  fclose(fs);
}


// Calculate total and ratio to test stats.
void PrintFilterAndPartialStats(bool print_found) {
  int total = 0;
  int total_to_test = 0;
  int partial_results = 0;
  int found = 4; // For d==0 (2,3,5,7)

  for (int d = START_DIGIT; d <= MAX_DIGITS; d++) {
    for (int a = 1; a <= 9; a++) {
      for (int b = 1; b <= 9; b += 2) {
        int status = is_prime[d][a][b];
        assert(status >= -2);

        if (b % 5 == 0 || (a + b) % 3 == 0) {
          continue;
        }

        total += 1;
        total_to_test += status == 0;
        partial_results += status == -2;
        found += status == -1;
      }
    }
  }

  int to_omit = 260; // Small easy values (less than 1e1000).
  if (print_found && found > to_omit) {
    printf("Showing %d partial results of %d found\n", found - to_omit, found);
    found = 4;
    for (int d = START_DIGIT; d <= MAX_DIGITS; d++) {
      for (int a = 1; a <= 9; a++) {
        for (int b = 1; b <= 9; b += 2) {
          if (is_prime[d][a][b] == -1) {
            found += 1;
            if (found >= to_omit) {
              cout << "\t" << found << " " << a << " * 10 ^ " << d << " + " << b << endl;
            }
          }
        }
      }
    }
  }

  printf("%d total, %d partial_results %d to test (%.3f remaining), (%.3f already tested)\n",
         total, partial_results, total_to_test, 1.0 * total_to_test / total, 1.0 * partial_results / (total_to_test + partial_results));
}


void WritePartialResult() {
  string file_name = "filter_" + to_string(START_DIGIT) + "_" + to_string(MAX_DIGITS) + ".partial";
  //cout << "\tWriting to: " << file_name << endl;
  ofstream fs(file_name, ios::app);

  // Write new_is_prime results to the file.
  int total_remaining = 0;
  int total_added = 0;
  for (int d = START_DIGIT; d <= MAX_DIGITS; d++) {
    for (int a = 1; a <= 9; a++) {
      for (int b = 1; b <= 9; b += 2) {
        if (is_prime[d][a][b] == 0) {
          total_remaining += 1;
        }

        if (new_status[d][a][b]) {
          assert(is_prime[d][a][b] <= -1);
          new_status[d][a][b] = false;
          total_added += 1;
          fs << d << ", " << a << ", " << b << ": " << is_prime[d][a][b] << endl;
        }
      }
    }
  }

  printf("%d total_remaining, added %d\n", total_remaining, total_added);
  fs.close();
}

float PredictCost() {
  float total_cost = 0;

  for (int d = START_DIGIT; d <= MAX_DIGITS; d++) {
    for (int a = 1; a <= 9; a++) {
      for (int b = 1; b <= 9; b++) {
        int status = is_prime[d][a][b];
        assert(status >= -2);
        if (status == 0) {
          total_cost += GetCost(a, d, b);
        }
      }
    }
  }

  return total_cost;
}

int main(void) {
  vector<string> seq;

  LoadFilter();
  LoadPartialResults();
  PrintFilterAndPartialStats(true /* print_found */);


  float predicted_cost = 1.0 / 0.0; // positiy infinity
  if (MAX_DIGITS > 400) {
    predicted_cost = PredictCost();
    if (START_DIGIT > 1) {
      printf("\tO(%d) - O(%d) ~= %.2f\n", MAX_DIGITS, START_DIGIT, predicted_cost);
    } else {
      printf("\tO(%d) ~= %.2f\n", MAX_DIGITS, predicted_cost);
    }
    cout << endl;
  }
  auto T0 = chrono::high_resolution_clock::now();

  float cost_done = 0;

  constexpr SMALL_D = 200;

  // Deal with fast (and numerous initial values)
  for (int d = START_DIGIT; d <= SMALL_D; d++) {
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
  for (int d = max(SMALL_D + 1, START_DIGIT); d <= MAX_DIGITS; d++) {
    cost_done += TestD(d);

    if (cost_done >= cost_done_print) {
      WritePartialResult();

      auto T1 = chrono::high_resolution_clock::now();
      auto elapsed_ms = chrono::duration_cast<chrono::milliseconds>(T1 - T0).count();
      float eta_ms = (predicted_cost / cost_done) * elapsed_ms;

      printf("Finished d: %d,  %.1f%% (%ld seconds, estimate: %.0f seconds)\n\n",
             d, 100 * cost_done / predicted_cost, elapsed_ms / 1000, eta_ms / 1000);
      cost_done_print += predicted_cost / status_prints;
    }
  }

  auto T1 = chrono::high_resolution_clock::now();
  auto primality_ms = chrono::duration_cast<chrono::milliseconds>(T1 - T0).count();
  cout << "\tCheck primality took " << primality_ms / 1000.0 << " seconds" << endl;
  cout << endl << endl;


  int count = 0;
  if (START_DIGIT == 1) {
    for (mpz_class t = 1; t < 10; t++) {
      if (mpz_probab_prime_p(t.get_mpz_t(), REPS)) {
        count += 1;
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

          // d = 1900 is around count == 290.
          if (count % 10 == 0 || d >= 1900) {
            cout << count << " " << a << " * 10^" << d << " + " << b << endl;
            WritePartialResult();
          }
        }
      }
    }
  }
}
