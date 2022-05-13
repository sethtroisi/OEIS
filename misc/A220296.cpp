#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <primesieve.hpp>
#include <string>
#include <vector>

using std::cout;
using std::endl;
using std::string;
using std::vector;

/**
 * Initially contains a single factor of f
 * In 2nd downward pass contains 0 if n is a n-digit home prime
 *   All primes are by default n-digit home primes
 *   re write all n-digit composite home primes to 0
 */
uint32_t *factor;

int main()
{
    /**
     * Needs 2GB for 9, 20GB for 10.
     * Could possible half using uint16_t for factor
     * Could reduce 1/3 or 2/3 by using wheel
     *
     * Could use bitset for is_home_prime
     * Then do factor in several passes
     */
    uint32_t max_d = 10;

    // Build factor table up to 10^d-1
    {
        uint64_t last = pow(10, max_d) - 1;
        uint64_t count = (last+1)/2;
        factor = (uint32_t*) malloc(sizeof(uint32_t) * count);
        std::fill(factor, factor + count, 0);

        primesieve::iterator it(2);
        uint64_t prime = it.next_prime();
        assert(prime == 3);
        for (; ; prime = it.next_prime()) {
            uint64_t p2 = prime * prime;
            if (p2 > last) break;

            // Only care about odd multiples
            for (uint64_t m = p2 / 2; m < count; m += prime) {
                factor[m] = prime;
            }
        }
        cout << "\tFactored up to prime: " << prime << endl;
    }

    vector<uint64_t> factors;
    for (uint64_t d = max_d; d >= 3; d--) {
        uint64_t first = pow(10, d-1);
        uint64_t last = pow(10, d) - 1;

        // Count of composite d-digit home primes
        uint64_t count_home_primes = 0;

        for (uint64_t lead_digit = 9; lead_digit >= 1; lead_digit--) {
            uint64_t first_with_lead = lead_digit * first;
            uint64_t last_with_lead = (lead_digit + 1) * first - 1;

            /**
             * hp iteration always produces a larger number
             * if you have a factor of 2 but don't start with a 1
             * that larger number will happen by having an extra digit
             */
            assert(last_with_lead % 2 == 1);
            uint64_t decr = lead_digit > 2 ? 2 : 1;

            // Not safe because we depend on factor from larger index being set for lower indexes.
            //#pragma omp parallel for schedule(dynamic, 128) private(factors) reduction(+:count_home_primes)
            for (uint64_t n = last_with_lead; n >= first_with_lead; n -= decr) {
                // Similiar logic to above, check for two factors of 2, 3
                if ((n & 3) == 0) continue;
                if (n % 9 == 0) continue;

                // Similiar logic to 2 above
                if (lead_digit > 3 && n % 3 == 0) continue;
                if (lead_digit > 5 && n % 5 == 0) continue;
                if (lead_digit > 7 && n % 7 == 0) continue;

                // If this is a prime continue
                if ((n & 1) && factor[n / 2] == 0)
                    continue;

                // for d=9 could use uint32_t and MAYBE for d=10 but not for d=11
                uint64_t t = n;
                factors.clear();

                // Normally this is a while loop, but n can only have one Factor.
                if ((t & 1) == 0) {
                    factors.push_back(2);
                    t /= 2;
                    //assert(t % 2 == 1); // See above
                }

                while (t > 1) {
                    uint64_t f = factor[t / 2];
                    if (f == 0) {
                        factors.push_back(t);
                        break;
                    } else {
                        factors.push_back(f);
                        t /= f;
                    }
                }

                if (factors.size() > d) {
                    continue;
                }

                std::sort(factors.begin(), factors.end());

                t = 1;
                uint64_t hp = 0;
                uint64_t pow_ten = 10;
                for (auto f : factors) {
                    t *= f;
                    while (f >= pow_ten) {
                        pow_ten *= 10;
                    }
                    hp *= pow_ten;
                    hp += f;
                }

                // cout << "\t" << n << ": " << hp << " |";
                // for (auto f : factors) cout << " " << f;
                // cout << endl;

                assert(t == n);
                assert(hp > n);

                if (hp <= last) {
                    /**
                     * hp will always be odd,
                     * Only even hp would be a powers of 2, where hp > last
                     */
                    assert(hp % 2 == 1);

                    // Check if hp continues to a d-digit-home-prime
                    if (factor[hp / 2] == 0) {
                        {
                            count_home_primes += 1;

                            // Don't have to store status of even numbers as they are never looked up
                            if (n % 2 == 1) {
                                //#pragma omp atomic write
                                factor[n / 2] = 0;
                            }
                        }
                    }
                }
            }
        }
        cout << d << " " << count_home_primes << endl;
    }
    return 0;
}
