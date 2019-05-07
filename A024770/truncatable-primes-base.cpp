#include <algorithm>
#include <gmp.h>
#include <gmpxx.h>
#include <iostream>
#include <vector>

#pragma omp declare reduction (merge : std::vector<mpz_class> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

using namespace std;

// Other bases result in absurd number of primes
vector<int> LEFT_BASES = {
2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34,
35, 37, 38, 39, 41, 43, 47, 49, 51, 53, 55, 59, 61, 65, 67,
71, 73, 79, 83, 89,
};

vector<mpz_class> SMALL_PRIMES = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37,
41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199,
211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293,
307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397,
401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491,
499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601,
607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701,
709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821,
823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929,
937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009
};

//#define LEFT false
#define LEFT true


long truncatable_primes(const int base) {
  long total = 0;

  vector<mpz_class> current = {};
  for (const auto p : SMALL_PRIMES) {
    if (p >= base) {
      break;
    }
    current.push_back(p);
  }


  mpz_class left_mult = 1;

  cout << "\t";
  for (int iteration = 1; !current.empty(); iteration += 1) {
    total += current.size();
    cout << "\t" << iteration << " " << current.size() << ", ";
    if (iteration % 5 == 0 || current.size() >= 1000000)
        cout << endl << "\t";

    left_mult *= base;

    vector<mpz_class> next;
    next.reserve(current.size());

    #pragma omp parallel for reduction(merge: next)
    for(auto it_v = current.begin(); it_v < current.end(); ++it_v) {

      #if LEFT
        for (int l = 1; l < base; l++) {
          mpz_class temp = l * left_mult + (*it_v);
      #else
        for (int r = 1; r < base; r += 2 - (base % 2)) {
          mpz_class temp = (*it_v) * base + r;
      #endif
          if (mpz_probab_prime_p(temp.get_mpz_t(), 25)) {
            next.push_back(temp);
          }
      }
    }

    swap(current, next);
  }
  cout << endl;

  return total;
};


int
main (void)
{
  for (int base = 2; base <= 100; base++) {
    #if LEFT
      if (find(LEFT_BASES.begin(), LEFT_BASES.end(), base) == LEFT_BASES.end())
        continue;
    #endif

    long result = truncatable_primes(base);
    cout << base << " " << result << endl;
  }
  return 0;
}
