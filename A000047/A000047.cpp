#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

#include <primesieve.hpp>
#include "../utils/count_special_primes.hpp"

using std::vector;

using std::cout;
using std::cerr;
using std::endl;


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
    const auto count_special_primes = get_special_prime_counts(
        n, r,
        /* start_prime= */ 3,
        [](uint64_t n) {
            uint64_t m = n % 8;
            return 2 * (n / 8) + (m >= 1) + (m >= 7);
        },
        [](uint64_t n) {
            uint64_t m = n % 8;
            return 2 * (n / 8) + (m >= 3) + (m >= 5);
        },
        [](uint64_t p) { uint8_t m = p & 7; return (m == 1) || (m == 7); }
    );

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
    count_in_ex = [&special_primes, &count_special_primes, &count_in_ex](uint64_t n, uint32_t pi) {
        if (n < special_primes[pi])
            return n;

        uint64_t count = n;

        // Handle p^3 < n
        for (; pi < special_primes.size(); pi++) {
            uint64_t p = special_primes[pi];
            uint64_t p2 = p * p;

            // This loop has been optimized see A000047.py, for clearer code
            uint64_t tn = n / p;
            if (p2 > tn)
                break;

            for (; ;) {
                if (tn < p) {
                    count -= tn;  // count_in_exp(tn, pi+1);
                    break;
                }
                count -= count_in_ex(tn, pi+1);

                // Have to add back all the counts of tn * p
                tn /= p;
                if (tn < p) {
                    count += tn;  // count_in_exp(tn, pi+1);
                    break;
                }
                count += count_in_ex(tn, pi+1);

                tn /= p;
            }
        }

        // Handle p^2 < n < p^3, only need to handle p^1 not p^3
        for (; pi < special_primes.size(); pi++) {
            uint32_t p = special_primes[pi];

            uint64_t tn = n / p;
            if (p > tn)
                break;

            count -= count_in_ex(tn, pi+1);
            // Have to add back all the counts of tn*r

            tn /= p;
            assert(tn < p);
            count += tn;  // count_in_exp(tn, pi+1);
        }

        // Handle primes > sqrt(n)
        uint64_t start_p = special_primes[pi];
        uint32_t first_m = n / start_p;
        //assert(start_p * start_p > n);
        assert(first_m < start_p);

        uint64_t last = n / (first_m + 1);
        uint64_t count_last = count_special_primes.at(last);
        for (uint32_t m = first_m; m > 0; m--) {
            // Count of number of primes with n / p == m
            //   -> Primes in the interval (n / (m + 1), n / m]
            uint64_t first = last;
            last  = n / m;

            if (m == first_m) {
                assert(first < last);
                assert(first <= special_primes.back());
                assert(first < start_p && start_p <= last);
            }

            uint64_t count_first = 0;
            if (first < start_p) {
                assert(m == first_m);
                assert(first <= special_primes.back());
                //count_first = upper_bound(special_primes.begin(), special_primes.end(), start_p - 1) - special_primes.begin();
                //assert(count_first == pi);
                count_first = pi;
            } else {
                count_first = count_last;
            }

            count_last = count_special_primes.at(last);

            assert(count_last >= count_first);
            // Is there a simplier version of this update that only uses count_last (or count_first?)
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
