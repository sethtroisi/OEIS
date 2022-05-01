#include "count_special_primes.hpp"

#include <cassert>
#include <chrono>
#include <cstdint>
#include <functional>
#include <unordered_map>
#include <utility>
#include <vector>

#include <primesieve.hpp>

using std::pair;
using std::vector;

/**
 * Get number of primes <= i for important values of i.
 * Assumes primes are in two congruence classes.
 *
 * See older code in A000047/A000205 for concrete examples
 * (number of primes % 8 in (5,7)) <= i for important values of i
 *
 * Adapted from Lucy_Hedgehog's post in Problem 10
 * https://projecteuler.net/thread=10;page=5#111677
 * https://math.stackexchange.com/a/2283829/87805
 */
Map<uint64_t, uint64_t>
get_special_prime_counts(
        uint64_t n, uint32_t r,
        uint32_t start_prime,
        std::function< uint64_t(uint64_t)> init_count_a,
        std::function< uint64_t(uint64_t)> init_count_b,
        std::function< bool(uint64_t)> is_group_a
) {
    auto start = std::chrono::high_resolution_clock::now();

    // Pair of how many numbers <= i of {form_a, form_b}
    // for i = n / 1, n / 2, ... n / r, n / r - 1, n / r - 2, ... 3, 2 , 1
    vector<pair<uint64_t, pair<uint64_t, uint64_t>>> counts_backing;
    {
        size_t size = r + n/r - 1;
        counts_backing.reserve(size);

        for(uint64_t i = 1; i <= r; i++) {
            uint64_t v = n / i;
            uint64_t c_a = init_count_a(v);
            uint64_t c_b = init_count_b(v);
            assert(c_a + c_b <= v);
            counts_backing.push_back({v, {c_a, c_b}});
        }

        for(uint32_t v = n / r - 1 ; v > 0; v--) {
            uint64_t c_a = init_count_a(v);
            uint64_t c_b = init_count_b(v);
            assert(c_a + c_b <= v);
            counts_backing.push_back({v, {c_a, c_b}});
        }
    }

    // Do calculation | 98% of the work is here
    {
        Map<uint64_t, pair<uint64_t, uint64_t>*> counts;
        counts.reserve(counts_backing.size() / 0.7);
        for (auto& [i, backing] : counts_backing) {
            counts[i] = &backing;
        }

        primesieve::iterator it(/* start= */ start_prime-1);
        uint64_t prime = it.next_prime();
        assert(prime == start_prime);
        for (; prime <= r; prime = it.next_prime()) {
            uint64_t p2 = prime * prime;

            auto [c_a, c_b] = *counts[prime-1];  // count of primes

            if (is_group_a(prime)) {
                for (auto& [v, u] : counts_backing) {
                    if (v < p2) break;

                    const auto& temp = counts[v / prime];
                    u.first  -= temp->first  - c_a;
                    u.second -= temp->second - c_b;
                }
            } else {
                for (auto& [v, u] : counts_backing) {
                    if (v < p2) break;

                    const auto& temp = counts[v / prime];
                    u.first  -= temp->second  - c_b;
                    u.second -= temp->first   - c_a;
                }
            }
        }
    }

    // Grab result in format we want
    Map<uint64_t, uint64_t> count_primes;
    count_primes.reserve(counts_backing.size() / 0.7);
    for (auto& [i, backing] : counts_backing) {
        //fprintf(stderr, "\t%lu\t%lu  %lu\n", i, backing.first, backing.second);
        count_primes[i] = backing.second;
    }

    {
        auto end = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double>(end - start).count();
        fprintf(stderr, "\tcount_special_primes(%lu) = %lu  (%.1f)\n",
                n, count_primes.at(n), elapsed);
    }
    return count_primes;
}
