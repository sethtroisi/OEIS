#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <utility>
#include <vector>

#include <primesieve.hpp>
#include "flat_hash_map.hpp"

using std::upper_bound;

using std::pair;
using std::vector;

using std::cout;
using std::cerr;
using std::endl;


template <class Key, class Val>
using Map = ska::flat_hash_map<Key, Val>;
//using Map = std::unordered_map<Key, Val>;

/**
 * Get number of primes % 8 == {3, 5} <= i for important values of i
 *
 * Adapted from Lucy_Hedgehog's post in Problem 10
 * https://projecteuler.net/thread=10;page=5#111677
 * https://math.stackexchange.com/a/2283829/87805
 */
Map<uint64_t, uint64_t>
get_three_five_prime_counts(uint64_t n, uint32_t r) {
    // Pair of how many numbers <= i of the form
    //  { 8*j + {1, 7}, 8*j + {3,5} }
    //    ^^^^^^^^^^^^ will includes the pseudoprime "1"
    vector<pair<uint64_t, pair<uint64_t, uint64_t>>> counts_backing;
    {
        // Convience vector so I don't have to duplicate c_a, c_b logic
        vector<uint64_t> V;
        {
            size_t size = r + n/r - 1;
            V.reserve(size);
            for(uint64_t i = 1; i <= r; i++) {
                V.push_back(n / i);
            }
            for(int32_t v = V[r-1] - 1; v > 0; v--) {
                V.push_back(v);
            }
            assert(V[0] == n);
            assert(V[size-1] == 1);
        }

        counts_backing.reserve(V.size());

        for (uint64_t i : V) {
            __uint128_t base = 2 * (i/8);
            char mod = i % 8;
            uint64_t c_a = base + (mod >= 1) + (mod >= 7);
            uint64_t c_b = base + (mod >= 3) + (mod >= 5);
            counts_backing.push_back({i, {c_a, c_b}});
        }
    }

    // Do calculation
    {
        Map<uint64_t, pair<uint64_t, uint64_t>* > counts;
        counts.reserve(counts_backing.size() / 0.7);
        for (auto& [i, backing] : counts_backing) {
            counts[i] = &backing;
        }


        primesieve::iterator it;
        uint64_t prime = it.next_prime();
        assert(prime == 2);
        // Only look at odd primes
        for (prime = it.next_prime(); prime <= r; prime = it.next_prime()) {
            uint64_t p2 = prime * prime;

            auto [c_a, c_b] = *counts[prime-1];  // count of primes: (8*k + {1,7}, 8*k + {3,5})

            if ((prime % 8) == 1 || (prime % 8 == 7)) {
                for (auto& [v, u] : counts_backing) {
                    if (v < p2) break;

                    pair<uint64_t, uint64_t> temp = *counts[v / prime];
                    u.first  -= temp.first  - c_a;
                    u.second -= temp.second - c_b;
                }
            } else {
                for (auto& [v, u] : counts_backing) {
                    if (v < p2) break;

                    pair<uint64_t, uint64_t> temp = *counts[v / prime];
                    u.first  -= temp.second  - c_b;
                    u.second -= temp.first   - c_a;
                }
            }
            /*
            if ((prime % 8) == 1 || (prime % 8 == 7)) {
                for (auto v : V) {
                    if (v < p2) break;

                    auto temp = counts[v / prime];
                    auto& u = counts[v];
                    u.first  -= temp.first  - c_a;
                    u.second -= temp.second - c_b;
                }
            } else {
                for (auto v : V) {
                    if (v < p2) break;

                    auto temp = counts[v / prime];
                    auto& u = counts[v];
                    u.first  -= temp.second  - c_b;
                    u.second -= temp.first   - c_a;
                }
            }
            */
        }
    }

    // Grab result in format we want
    Map<uint64_t, uint64_t> count_primes;
    count_primes.reserve(counts_backing.size() / 0.7);
    for (auto& [i, backing] : counts_backing) {
        count_primes[i] = backing.second;
    }
    return count_primes;
}


uint64_t A000047_final(size_t bits) {
    uint64_t n = 1ul << bits;
    uint64_t r = sqrt(n) + 1;
    while (r*r > n) {
        r--;
    }
    assert(r*r <= n);
    assert((r+1) * (r+1) > n);
    assert(r < std::numeric_limits<uint32_t>::max());


    // 10-50% of time is building special prime counts.
    const auto count_special_primes = get_three_five_prime_counts(n, r);
    cerr << "\tcount_special_primes(2^" << bits << ") = " << count_special_primes.at(n) << endl;
    // return count_special_primes.at(n);

    // Build list of special primes p % 8 == {3, 5}
    // Only interested in these primes to odd powers

    vector<uint32_t> special_primes;
    {
        size_t past = 0;
        primesieve::iterator it;
        for (uint64_t prime = it.next_prime(); past < 2; prime = it.next_prime()) {
            if (prime % 8 == 3 || prime % 8 == 5) {
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

    std::function<uint64_t(uint64_t, uint64_t)> count_in_ex;
    count_in_ex = [&r, &special_primes, &count_special_primes, &count_in_ex](uint64_t n, uint32_t pi) {
        if (n < special_primes[pi])
            return n;

        uint64_t count = n;
        for (; pi < special_primes.size(); pi++) {
            uint64_t p = special_primes[pi];
            uint64_t p2 = p * p;
            if (p2 > n)
                break;

            // This loop has been optimized see A000047.py, for clearer code
            uint64_t tn = n / p;
            for (; ;) {
                if (tn < p) {
                    count -= tn;  // count_in_exp(tn, pi+1);
                    break;
                }
                count -= count_in_ex(tn, pi+1);

                // Have to add back all the counts of tn*r
                tn /= p;
                if (tn < p) {
                    count += tn;  // count_in_exp(tn / p, pi+1);
                    break;
                }
                count += count_in_ex(tn, pi+1);

                tn /= p;
            }
        }

        // Handle primes > sqrt(n)
        uint64_t start_p = (pi < special_primes.size()) ? special_primes[pi] : (special_primes.back() + 1);
        assert(start_p * start_p > n);
        uint32_t first_m = n / start_p;

        for (uint32_t m = first_m; m > 0; m--) {
            // Count of number of primes with n / p == m
            //   -> Primes in the interval (n / (m + 1), n / m]

            uint64_t first = n / (m + 1);
            uint64_t last  = n / m;

            if (m == first_m) {
                assert(first < last);
                assert(first <= special_primes.back());
                assert(first < start_p <= last);
            }

            uint64_t count_first = 0;
            if (first < start_p) {
                assert(m == first_m);
                assert(first <= special_primes.back());
                //count_first = upper_bound(special_primes.begin(), special_primes.end(), start_p - 1) - special_primes.begin();
                //assert(count_first == pi);
                count_first = pi;
            } else {
                count_first = count_special_primes.at(first);
            }

            uint64_t count_last = count_special_primes.at(last);

            assert(count_last >= count_first);
            count -= m * (count_last - count_first);
        }

        return count;
    };

    return count_in_ex(n, 0);
}

int main(int argc, char** argv) {
    size_t bits = 30;
    if (argc == 2) {
        bits = atoi(argv[1]);
    }

    uint64_t count = A000047_final(bits);
    cout << "A000047(" << bits << ") = " << count << endl;
    return 0;
}
