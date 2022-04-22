#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <unordered_map>
#include <utility>
#include <vector>

#include <primesieve.hpp>

using std::upper_bound;

using std::pair;
using std::unordered_map;
using std::vector;

using std::cout;
using std::cerr;
using std::endl;


/**
 * Get number of primes % 8 == {3, 5} <= i for important values of i
 *
 * Adapted from Lucy_Hedgehog's post in Problem 10
 * https://projecteuler.net/thread=10;page=5#111677
 * https://math.stackexchange.com/a/2283829/87805
 */
unordered_map<uint64_t, uint64_t>
get_three_five_prime_counts(uint64_t n, uint32_t r) {
    vector<uint64_t> V;
    {
        size_t size = r + n/r - 1;
        V.reserve(size);
        for(uint64_t i = 1; i <= r; i++) {
            V.push_back(n / i);
        }
        // Walk downwards, moving forwards
        for(int32_t v = V[r-1] - 1; v > 0; v--) {
            V.push_back(v);
        }
        assert(V[0] == n);
        assert(V[size-1] == 1);
    }

    // Pair of how many numbers <= i of the form
    //  { 8*j + {1, 7}, 8*j + {3,5} }
    unordered_map<uint64_t, pair<uint64_t, uint64_t>> counts;
    counts.reserve(V.size());
    for (uint64_t i : V) {
        uint64_t base = 2 * (i/8);
        char mod = i % 8;
        counts[i] = {
            base + (mod >= 1) + (mod >= 7),
            base + (mod >= 3) + (mod >= 5)
        };
    }

    primesieve::iterator it;
    uint64_t prime = it.next_prime();
    assert(prime == 2);
    // Only look at odd primes
    for (prime = it.next_prime(); prime <= r; prime = it.next_prime()) {
        uint64_t p2 = prime * prime;

        //auto [c_a, c_b] = counts[prime-1];  // count of primes: (8*k + {1,7}, 8*k + {3,5})
        uint64_t c_a = counts[prime-1].first;
        uint64_t c_b = counts[prime-1].second;

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
    }

    // c_a also includes the pseudoprime "1"

    unordered_map<uint64_t, uint64_t> count_primes;
    for (const auto& it : counts) {
        count_primes[it.first] = it.second.second;
    }
    return count_primes;
}


uint64_t A000047_final(size_t bits) {
    uint64_t n = 1ul << bits;
    uint32_t r = sqrt(n) + 1;
    while (r*r > n) {
        r--;
    }
    assert(r*r <= n);

    // Need slightly more than sqrt(r) primes
    // primes = get_prime_array(r + 100)
    // print(f"Primes({len(primes)}) {primes[0]} to {primes[-1]}")

    // 40-70% of time is this call
    const auto count_special_primes = get_three_five_prime_counts(n, r);
    cerr << "\tcount_special_primes(2^" << bits << ") = " << count_special_primes.at(n) << endl;
    // return count_special_primes.at(n);

    // Build list of special primes p % 8 == {3, 5}
    // Only interested in these primes to odd powers

    vector<uint32_t> special_primes;
    primesieve::iterator it;
    for (uint32_t prime = it.next_prime(); prime <= r || special_primes.back() < r; prime = it.next_prime()) {
        if (prime % 8 == 3 || prime % 8 == 5)
            special_primes.push_back(prime);
    }
    assert(special_primes.back() > r);  // Need one past r

    std::function<uint64_t(uint64_t, uint64_t)> count_in_ex;
    count_in_ex = [&r, &special_primes, &count_special_primes, &count_in_ex](uint64_t n, uint32_t pi) {
        if (n < special_primes[pi])
            return n;

        uint64_t count = n;
        for (; pi < special_primes.size(); pi++) {
            uint32_t p = special_primes[pi];
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
        assert(pi < special_primes.size());
        uint64_t start_p = special_primes[pi];
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

            uint64_t count_first;
            if (first < start_p) {
                assert(m == first_m);
                assert(first <= special_primes.back());
                count_first = upper_bound(special_primes.begin(), special_primes.end(), start_p - 1) - special_primes.begin();
            } else {
                count_first = count_special_primes.at(first);
                // Nice double check of special_prime code
                //uint64_t test = upper_bound(special_primes.begin(), special_primes.end(), first) - special_primes.begin();
                //assert(count_first == test);
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
