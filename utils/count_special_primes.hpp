#include <cstdint>
#include <functional>
#include <unordered_map>
#include <vector>

template <class Key, class Val>
//#include "flat_hash_map.hpp"
//using Map = ska::flat_hash_map<Key, Val>;
using Map = std::unordered_map<Key, Val>;

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
);