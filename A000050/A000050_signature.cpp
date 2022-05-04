#include <algorithm>
#include <cassert>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <vector>

#include <primesieve.hpp>
#include "../utils/count_special_primes.hpp"


using std::vector;

using std::cout;
using std::cerr;
using std::endl;


uint64_t A000050_final(size_t bits) {
    uint64_t n = 1ul << bits;
    uint64_t r = sqrt(n) + 1;
    while (r*r > n) {
        r--;
    }
    assert(r*r <= n);
    assert((r+1) * (r+1) > n);
    assert(r < std::numeric_limits<uint32_t>::max());

    auto start = std::chrono::high_resolution_clock::now();

    // 30-80% of time is building special prime counts.
    const auto count_special_primes = get_special_prime_counts_vector(
        n, r,
        /* start_prime= */ 3,
        [](uint64_t n) { return (n / 4) + ((n % 4) >= 1); },
        [](uint64_t n) { return (n / 4) + ((n % 4) >= 3); },
        [](uint64_t p) { return (p & 3) == 1; }
    );

    const auto special_counts = count_special_primes.size();
    assert(special_counts == (n/r-1) + r);

    std::function<uint32_t(uint64_t)> count_special_index
      = [special_counts, n, r](uint64_t v) {
      uint64_t i = (v <= r) ? v - 1 : (special_counts - n / v);
      return i;
    };

    // Build list of special primes p % 4 == 3}
    // Only interested in these primes to odd powers
    vector<uint32_t> special_primes;
    {
        size_t past = 0;
        primesieve::iterator it;
        for (uint64_t prime = it.next_prime(); past < 2; prime = it.next_prime()) {
            if (prime % 4 == 3) {
                special_primes.push_back(prime);
                past += prime > r;
            }
        }
        assert(special_primes[special_primes.size() - 2] > r);  // Need two past r
        cerr << "\tPrimes(" << special_primes.size() << ") = "
            << special_primes[0] << " ... "
            << special_primes[special_primes.size() - 3] << ", "
            << special_primes[special_primes.size() - 2] << ", "
            << special_primes[special_primes.size() - 1] << endl;
    }

    std::function<uint64_t(uint64_t, uint32_t)> count_in_ex_large
      = [&special_primes, &count_special_primes, &count_special_index](uint64_t n, uint32_t pi) {
        // Generally happens when primes[pi-1] <= n < primes[pi]
        if (n < special_primes[pi]) {
            return n;
        }

        uint64_t count = n;

        // Handle p <= n < p^2
        uint64_t start_p = special_primes[pi];
        assert(start_p * start_p > n);

        uint32_t first_m = n / start_p;
        assert(first_m < start_p);

        uint64_t last = n / (first_m + 1);
        uint64_t count_last = 0;

        // Determine start prime of first loop of sqrt code.
        if (last < start_p) {
          assert(last <= special_primes.back());
          count_last = pi;
        } else {
          //count_last = count_special_primes.at(last);
          count_last = count_special_primes[count_special_index(last)];
        }

        for (uint32_t m = first_m; m > 0; m--) {
            // Count of number of primes with n / p == m
            //   -> Primes in the interval (n / (m + 1), n / m]
            uint64_t first = last;
            last  = n / m;
            uint64_t count_first = count_last;

            //count_last = count_special_primes.at(last);
            count_last = count_special_primes[count_special_index(last)];

            assert(count_last >= count_first);
            // Is there a simplier version of this update that only uses count_last (or count_first?)
            count -= m * (count_last - count_first);
        }
        return count;
    };

    std::function<uint64_t(uint64_t, uint32_t, uint8_t)> count_in_ex;
    count_in_ex = [&special_primes, &count_in_ex, &count_in_ex_large]
                      (uint64_t n, uint32_t pi, uint8_t root) {
        if (n < special_primes[pi])
            return n;

        uint64_t count = 0;

        if (root >= 4) {
          // Handle p where p^4 < n
          for (; pi < special_primes.size(); pi++) {
              uint64_t p = special_primes[pi];
              uint64_t p2 = p * p;

              // Do I need to worry about overflow (on first index) for recursive calls?
              // Could also trigger if overflow happens
              if (p2 > n / p2)
                  break;

              // This loop has been optimized see A000047.py, for clearer code
              uint64_t tn = n / p;

              for (; ;) {
                  if (tn < p) {
                      count -= tn;  // count_in_exp(tn, pi+1);
                      break;
                  }
                  count -= count_in_ex(tn, pi+1, root);

                  // Have to add back all the counts of tn * p
                  tn /= p;
                  if (tn < p) {
                      count += tn;  // count_in_exp(tn, pi+1);
                      break;
                  }
                  count += count_in_ex(tn, pi+1, root);

                  tn /= p;
              }
          }
        }

        if (root >= 3) {
          // Handle p where p^3 <= n < p^4
          // p^1 has recursion twice
          // p^2 has trivial recursion
          // p^3 has no recursion
          // Handle p^2 < n < p^3, only need to handle p^1 not p^3
          for (; pi < special_primes.size(); pi++) {
              uint32_t p = special_primes[pi];
              uint64_t p2 = p * p;

              uint64_t tn = n / p;
              if (p2 > tn)
                  break;

              // p^3 < n < p^4  ->  p^2 <= tn < p^3

              // TODO could possible pass a "skip first loop flag" (or call a
              // method that handles large p)
              count -= count_in_ex(tn, pi+1, 2);

              tn /= p; // p <= tn < p^2
              // For each of these primes we are going to run the final loop logic
              count += count_in_ex_large(tn, pi+1);

              tn /= p; // 1 <= tn < p
              assert(tn < p);
              count -= tn;
          }
        }
//*/

        if (root >= 2) {
          // Handle p^2 <= n < p^3, only need to handle p^1 not p^3
          for (; pi < special_primes.size(); pi++) {
              uint32_t p = special_primes[pi];

              uint64_t tn = n / p;
              if (p > tn)
                  break;

              count -= count_in_ex_large(tn, pi+1);
              // Have to add back all the counts of tn*r

              tn /= p;
              assert(tn < p);
              count += tn;  // count_in_exp(tn, pi+1);
          }
        }

        // Handles adding n to count.
        count += count_in_ex_large(n, pi);

        return count;
    };

    uint64_t count = count_in_ex(n, 0, 100);

    {
        auto end = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double>(end - start).count();
        printf("| %2lu | %-13lu | %-14lu |              |          | %-7.2f |\n",
                bits, count, count_special_primes.back(), elapsed);
    }
    return count;
}

int main(int argc, char** argv) {
    size_t bits = 30;
    if (argc == 2) {
        bits = atoi(argv[1]);
    }

    uint64_t count = A000050_final(bits);
    cout << "A000050(" << bits << ") = " << count << endl;
    return 0;
}
