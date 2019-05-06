#include <gmp.h>
#include <gmpxx.h>
#include <iostream>
#include <vector>

#pragma omp declare reduction (merge : std::vector<mpz_class> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

using namespace std;

int
main (void)
{
  long total = 0;
  vector<mpz_class> current = {2, 3, 5, 7};
//  vector<mpz_class> current = {11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};

  mpz_class left_mult = current.back() >= 10 ? 10 : 1;

  for (int iteration = 1; !current.empty() && iteration <= 1000; iteration += 2) {
    total += current.size();
    cout << iteration << " " << total << " " << current.size() << endl;

    left_mult *= 100;

    vector<mpz_class> next;
    next.reserve(current.size());

    #pragma omp parallel for reduction(merge: next)
    for(auto it_v = current.begin(); it_v < current.end(); ++it_v) {
      for (int l = 1; l <= 9; l++) {
        mpz_class left_temp = l * left_mult + 10 * (*it_v);
        for (int r = 1; r <= 9; r += 2) {
          if (r == 5) continue;

          mpz_class temp = left_temp + r;
          if (mpz_probab_prime_p(temp.get_mpz_t(), 25)) {
            next.push_back(temp);
          }
        }
      }
    }

    swap(current, next);
  }

  return 0;
}
