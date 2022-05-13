#include <cassert>
#include <cstdint>
#include <iostream>
#include <primesieve.hpp>

int main()
{
    uint64_t divisor = 1;
    for (uint64_t n = 1; n < 15; n++) {
        divisor *= 7;
        uint64_t pow_ten = 10;
        uint64_t pow_ten_mod = 10;

        uint64_t smarandache_wellin_mod = 2;
        primesieve::iterator it(2);
        uint64_t prime;
        for (; smarandache_wellin_mod != 0;) {
            prime = it.next_prime();
            if (prime >= pow_ten) {
                pow_ten *= 10;
                pow_ten_mod = pow_ten % divisor;
                std::cout << "\t" << prime << " now " << pow_ten << " % divisor = " << pow_ten_mod << std::endl;
                __uint128_t temp;
                // Verify no overflow possible from mod * 10^digits + prime
                //assert(__builtin_umulll_overflow(pow_ten_mod, divisor, &temp) == 0);
                //assert(temp / pow_ten_mod == divisor);
                //assert(__builtin_umulll_overflow(temp, pow_ten, &temp) == 0);
                temp = (__uint128_t) pow_ten_mod * divisor;
                assert(temp / pow_ten_mod == divisor);
                assert((temp + 10 * pow_ten) > temp);
            }

            //smarandache_wellin_mod *= pow_ten_mod;
            //smarandache_wellin_mod += prime;
            //smarandache_wellin_mod %= divisor;

            __uint128_t temp = ((__uint128_t) smarandache_wellin_mod * pow_ten_mod + prime) % divisor;
            smarandache_wellin_mod = temp;
        }
        std::cout << n << " " << prime << std::endl;
    }
    return 0;
}
