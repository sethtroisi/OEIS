#include "count_special_primes.hpp"
#include <sys/types.h>

#include <cassert>
#include <chrono>
#include <cstdint>
#include <functional>
#include <utility>
#include <vector>

#include <primesieve.hpp>

using std::pair;
using std::vector;

/**
 * First attempt at parallelizing failed
 *
 * Let
 *   i3 = iroot<3>(n)
 *   i2 = r = iroot<2>(n)
 *
 * Easy idea is to only parallel when p > i3.
 *   * Nice because only top 1/3 of counts is updated
 *   * Never any interaction with counts[n] depending on counts[n/p] which hasn't been updated yet
 * BUT
 *   * only a small amount of time is processing these large primes (think about triangle numbers)
 *   * counts[n] gets cobbered by each seperate p trying to update it
 *     * could work around this by inverting the two loops
 *
 * Not quite sure what's the best route from here
 */


vector<pair<uint64_t, pair<uint64_t, uint64_t>>>
__get_special_prime_counts(
        uint64_t n, uint32_t r,
        uint32_t start_prime,
        std::function< uint64_t(uint64_t)> init_count_a,
        std::function< uint64_t(uint64_t)> init_count_b,
        std::function< bool(uint64_t)> is_group_a
) {
    // Pair of how many numbers <= i of {form_a, form_b}
    // for i = 1, 2, ..., n/r, n/(r-1), n/(r-2), ... n/3 n/2 n/1

    // NOTE: Can technically drop "v" from vector (requires compute v in inner loop)
    vector<pair<uint64_t, pair<uint64_t, uint64_t>>> counts_backing;
    {
        size_t size = r + n/r - 1;
        counts_backing.reserve(size);
        // 1, 2, ... n / r - 1
        for(uint32_t v = 1; v < (n / r); v++) {
            uint64_t c_a = init_count_a(v);
            uint64_t c_b = init_count_b(v);
            assert(c_a + c_b <= v);
            counts_backing.push_back({v, {c_a, c_b}});
        }

        // n/r, n/(r-1), n/(r-2), ... n/3 n/2 n/1
        for(uint64_t i = r; i >= 1; i--) {
            uint64_t v = n / i;
            uint64_t c_a = init_count_a(v);
            uint64_t c_b = init_count_b(v);
            assert(c_a + c_b <= v);
            counts_backing.push_back({v, {c_a, c_b}});
        }
    }

    const auto length = ((n/r)-1) + r;
    assert(counts_backing.size() == length);

    // Do calculation | 98% of the work is here
    {
        primesieve::iterator it(/* start= */ start_prime-1);
        uint64_t prime = it.next_prime();
        assert(prime == start_prime);
        for (; prime <= r; prime = it.next_prime()) {
            uint64_t p2 = prime * prime;

            auto& [v, c__] = counts_backing[prime-2];
            assert(v == prime-1);
            auto [c_a, c_b] = c__;

            // Index of last term < p2
            uint64_t stop_i = ((p2-1) < r) ? (p2-2) : (length - ((n-1) / p2 + 1));
            assert(counts_backing[stop_i].first < p2);
            assert(counts_backing[stop_i+1].first >= p2);

            if (is_group_a(prime)) {
                for (size_t i = counts_backing.size() - 1; i > stop_i; i--) {
                    auto& [v, u] = counts_backing[i];
                    assert(v >= p2);

                    // Only need v for this t
                    uint64_t t = v / prime;
                    // t = (n / i) / prime
                    size_t index = (t < r) ? (t-1) : (length - (n / t));
                    assert(counts_backing[index].first == t);

                    const auto temp = counts_backing[index].second;
                    u.first  -= temp.first  - c_a;
                    u.second -= temp.second - c_b;
                }
            } else {
                for (size_t i = counts_backing.size() - 1; i > stop_i; i--) {
                    auto& [v, u] = counts_backing[i];
                    assert(v >= p2);

                    uint64_t t = v / prime;
                    size_t index = (t < r) ? (t-1) : (length - (n / t));
                    assert(counts_backing[index].first == t);

                    const auto temp = counts_backing[index].second;
                    u.first  -= temp.second  - c_b;
                    u.second -= temp.first   - c_a;
                }
            }
        }
    }
    return counts_backing;
}


// /**
//  * Get number of primes <= i for important values of i.
//  * Assumes primes are in two congruence classes.
//  *
//  * See older code in A000047/A000205 for concrete examples
//  * (number of primes % 8 in (5,7)) <= i for important values of i
//  *
//  * Adapted from Lucy_Hedgehog's post in Problem 10
//  * https://projecteuler.net/thread=10;page=5#111677
//  * https://math.stackexchange.com/a/2283829/87805
//  */
// Map<uint64_t, uint64_t>
// get_special_prime_counts(
//         uint64_t n, uint32_t r,
//         uint32_t start_prime,
//         std::function< uint64_t(uint64_t)> init_count_a,
//         std::function< uint64_t(uint64_t)> init_count_b,
//         std::function< bool(uint64_t)> is_group_a
// ) {
//     auto start = std::chrono::high_resolution_clock::now();
//
//     auto counts_backing = __get_special_prime_counts(
//         n, r, start_prime,
//         init_count_a, init_count_b, is_group_a);
//
//     // Grab result in format we want
//     Map<uint64_t, uint64_t> count_primes;
//     count_primes.reserve(counts_backing.size() / 0.7);
//     for (auto& [i, backing] : counts_backing) {
//         //fprintf(stderr, "\t%lu\t%lu  %lu\n", i, backing.first, backing.second);
//         count_primes[i] = backing.second;
//     }
//
//     {
//         auto end = std::chrono::high_resolution_clock::now();
//         double elapsed = std::chrono::duration<double>(end - start).count();
//         fprintf(stderr, "\tcount_special_primes(%lu) = %lu  (%.2f)\n",
//                 n, count_primes.at(n), elapsed);
//     }
//     return count_primes;
// }


/**
 * Get number of primes <= i for important values of i.
 * Returns result in vector
 *
 * See get_special_prime_counts_map
 */
vector<uint64_t>
get_special_prime_counts_vector(
        uint64_t n, uint32_t r,
        uint32_t start_prime,
        std::function< uint64_t(uint64_t)> init_count_a,
        std::function< uint64_t(uint64_t)> init_count_b,
        std::function< bool(uint64_t)> is_group_a
) {
    auto start = std::chrono::high_resolution_clock::now();

    auto counts_backing = __get_special_prime_counts(
    //auto counts_backing = __get_special_prime_counts(
        n, r, start_prime,
        init_count_a, init_count_b, is_group_a);

    assert(counts_backing[0].first == 1);      // First value should be 1
    assert(counts_backing.back().first == n);  // Last value should be n

    /**
     * Reduces 3x uint64 to 1x uint64
     * BUT temporarily increases memory pressure
     *
     * Possibly I can read/erase from back, then reserve
     */
    vector<uint64_t> count_primes;
    count_primes.resize(counts_backing.size());
    for (size_t i = 0; i < counts_backing.size(); i++) {
        count_primes[i] = counts_backing[i].second.second;
    }
    // counts should be strictly increasing
    assert(is_sorted(count_primes.begin(), count_primes.end()));

    {
        auto end = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double>(end - start).count();
        fprintf(stderr, "\tcount_special_primes(%lu) = %lu  (%.2f)\n",
                n, count_primes.back(), elapsed);
    }
    return count_primes;
}


/**
 * Count population of 2^n of quadratic form generated by [2][p][q]^2
 *
 * TODO: can I calculate this directly using a variation of
 * get_special_prime_counts. Possibly computing "number less than n without
 * a duplicate factor".
 */
uint64_t count_population_quadratic_form(
        size_t bits,
        uint32_t start_prime,
        uint32_t add_to_special_primes,
        std::function< uint64_t(uint64_t)> init_count_a,
        std::function< uint64_t(uint64_t)> init_count_b,
        std::function< bool(uint64_t)> is_group_a
) {
    uint64_t n = 1ul << bits;
    uint64_t r = sqrt(n) + 1;
    while (r*r > n) {
        r--;
    }
    assert(r*r <= n);
    assert((r+1) * (r+1) > n);
    assert(r < std::numeric_limits<uint32_t>::max());

    auto start = std::chrono::high_resolution_clock::now();

    // 40-60% of time is building special prime counts.
    const auto count_special_primes = get_special_prime_counts_vector(
        n, r,
        start_prime,
        init_count_a,
        init_count_b,
        is_group_a);

    const auto special_counts = count_special_primes.size();
    assert(special_counts == (n/r-1) + r);

    // TODO does (v <= r) work when n does not evenly divide r
    // TODO try inlining
    std::function<uint32_t(uint64_t)> count_special_index
      = [special_counts, n, r](uint64_t v) {
      return (v <= r) ? v - 1 : (special_counts - n / v);
    };

    // Only interested in these primes to odd powers
    vector<uint32_t> special_primes;
    {
        if (add_to_special_primes) {
          special_primes.push_back(add_to_special_primes);
        }

        size_t past = 0;
        primesieve::iterator it(/* start= */ start_prime-1);
        for (uint64_t prime = it.next_prime(); past < 2; prime = it.next_prime()) {
            if (!is_group_a(prime)) {
                special_primes.push_back(prime);
                past += prime > r;
            }
        }
        assert(special_primes[special_primes.size() - 2] > r);  // Need two past r
        fprintf(stderr, "\tPrimes(%lu) = %u %u ... %u, %u, %u\n",
            special_primes.size(),
            special_primes[0],
            special_primes[1],
            special_primes[special_primes.size() - 3],
            special_primes[special_primes.size() - 2],
            special_primes[special_primes.size() - 1]);
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
          count_last = count_special_primes[count_special_index(last)];
        }

        for (uint32_t m = first_m; m > 0; m--) {
            // Count of number of primes with n / p == m
            //   -> Primes in the interval (n / (m + 1), n / m]
            // uint64_t first = last;
            last  = n / m;
            uint64_t count_first = count_last;

            count_last = count_special_primes[count_special_index(last)];

            assert(count_last >= count_first);
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
              assert(tn >= p);

              // For each of these primes run the final loop logic
              count += count_in_ex_large(tn, pi+1);

              tn /= p; // 1 <= tn < p
              assert(tn < p);
              count -= tn;
          }
        }

        if (root >= 2) {
          // Handle p^2 <= n < p^3, only need to handle p^1 not p^3
          for (; pi < special_primes.size(); pi++) {
              uint32_t p = special_primes[pi];

              uint64_t tn = n / p;
              if (p > tn)
                  break;

              assert(tn >= p);
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
        printf("| %2lu | %-18lu | %-18lu | %-8.2f |\n",
                bits, count, count_special_primes.back(), elapsed);
    }
    return count;
}
