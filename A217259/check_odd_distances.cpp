// g++ -g -O3 --std=c++17 -Werror -Wall check_odd_distances.cpp

#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <vector>

using std::vector;


vector<uint32_t> gen_primes(uint32_t min, uint32_t n) {
    vector<bool> status(n+1, 1);

    for (uint32_t p = 3; p*p <= n; p += 2) {
        if (status[p]) {
            for (uint32_t m = p*p; m <= n; m += p) {
                status[m] = 0;
            }
        }
    }

    vector<uint32_t> primes;
    if (2 >= min) {
        primes.push_back(2);
    }
    for (uint32_t p = (min | 1); p <= n; p += 2) {
        if (status[p]) {
            if (p >= min)
                primes.push_back(p);
        }
    }
    return primes;
}


// Only works for uint32_t inputs
uint64_t powMod(uint64_t base, uint64_t exp, uint64_t modulus)
{
    base = base % modulus;

    uint64_t result = 1;
    while (exp > 0) {
        if (exp & 1) {
            result = (result * base) % modulus;
        }
        base = (base * base) % modulus;
        exp >>= 1;
    }
    return result;
}

uint64_t tonelli_shanks(uint64_t p, uint64_t n) {
    // find R such that R*R mod p = n
    //   (or if no R exists return None)

    assert( p < 1'000'000'000 );
    if (n < 0)
        n += p;

    if (n == 0)
        return 0;

    int32_t legrendre = powMod(n, (p-1) / 2, p);
    assert (legrendre == 1);

    uint64_t Q = p - 1;
    uint64_t S = 0;
    while ((Q & 1) == 0) { // Q % 2 == 0
        Q >>= 1; // Q /= 2
        S += 1;
    }

    if (S == 1) {
        // p % 3 == 4
        return powMod(n, (p+1) / 4, p);
    }

    uint64_t z = 2;
    while (powMod(z, (p-1) / 2, p) != (p-1))
        z += 1;
    // z is now an element with order p-1

    //cout << "Q: " << Q << " " << S << " " << z << endl;

    uint64_t c = powMod(z, Q, p);
    uint64_t R = powMod(n, (Q+1)/2, p);
    uint64_t t = powMod(n, Q, p);
    uint64_t M = S;
    assert (M < 32);

    while (t != 1) {
        int i = 1;
        uint64_t tempT = (t*t) % p;
        while (tempT != 1) {
            i += 1;
            tempT = (tempT*tempT) % p;
        }
        // invariant is that pow(t, 2 ** i, p) == 1
        uint64_t b = powMod(c, 1 << (M-i-1), p);
        uint64_t bb = (b*b) % p;
        R = (R*b) % p;
        t = (t*bb) % p;
        c = bb;
        M = i;
    }

    return R;
}

typedef struct {
  uint64_t product;
  uint64_t sigma;
} num;


uint32_t isqrt(uint64_t n) {
    uint64_t t = sqrt(n) + 1;
    assert(t * t > n);
    while (t * t > n) { t--; }
    assert(t * t <= n);
    return t;
}



/**
 * Factor a^2 - offset  <=>    a^2 = offset mod p
 * offset is 0 or a small odd number and can be negative
 */
const vector<uint64_t> factor_offset(
        const uint32_t N,
        const bool twice_square,
        const int32_t offset,
        const vector<uint32_t> &primes) {
    assert((offset == 0) || ((offset & 1) == 1));
    uint32_t square_mult = twice_square ? 2 : 1;

    // [i'th term product of factors, i'th term partial sigma product, i+1'th term ..., i+1'th term, ...]
    vector<num> status(N+1, {1, 1});
    assert(status[5].product == 1 && status[5].sigma == 1);

    // Find first {2,1} * n^2 - offset > 2, sigma(i) below 2 doesn't make sense
    const uint32_t start_index = offset <= 0 ? 2 : (isqrt((2 + offset + square_mult - 1) / square_mult - 1) + 1);
    //printf("\t\tFactoring %u * n^2 %+d start at n=%u\n", 1 + twice_square, -offset, start_index);
    assert(((int32_t) square_mult * start_index * start_index - offset) >= 2);

    uint64_t factor_count = 0;

    /**
     * Remove prime (and possible prime^2, prime^3)
     * From status[base], status[base+prime], status[base+2*prime], ...
     */
    std::function remove_factor{
        [&](uint32_t base, uint32_t prime) {
            while (base < start_index) {
                base += prime;
            }
            //assert(((uint64_t) square_mult * base * base - offset) % prime == 0);
            //printf("\t\t%u at index %u\n", prime, base);

            for (uint32_t index = base; index <= N; index += prime) {
                uint64_t num = (uint64_t) square_mult * index * index - offset;
                //assert(num % prime == 0);
                num /= prime;

                uint64_t pp = prime;
                uint32_t exponent = 1;

                // TODO check if same 20-40% of time is spent checking for p^2 as in python
                while (num) {
                    uint64_t temp = num / prime;
                    if (temp * prime < num)
                        break;

                    num = temp;
                    pp *= prime;
                    exponent += 1;
                }

                factor_count += 1;
                //assert(save % pp == 0);
                status[index].product *= pp;
                status[index].sigma *= (prime * pp - 1) / (prime - 1);
                //printf("\t\t%u at index %u -> %lu | %lu & %lu\n", prime, base, save,
                //        status[index].product, status[index].sigma);
            }
        }
    };

    assert(primes[0] == 3);
    for (const auto prime : primes) {
        // offset can cause issues when negative.
        uint32_t residual = (prime + (offset % ((int32_t) prime))) % prime;
        if (twice_square) {
            residual = ((uint64_t) residual * (prime+1)/2) % prime;
        }

        if (residual == 0) {
            remove_factor(0, prime);
        } else {
            // Check if quadratic residue exists
            int32_t legendre = powMod(residual, (prime-1)/2, prime);

            if (legendre == 1) {
                // Find "square root" of offset, using tonelli-shanks
                uint32_t base = tonelli_shanks(prime, residual);
                assert(1 <= base && base < prime);

                uint32_t other = prime - base;
                assert(base != other);

                remove_factor(base, prime);
                remove_factor(other, prime);
            }
        }
    }

    //printf("\t%'7lu prime factors in %sn^2 %+d\n", factor_count, twice_square ? "2*" : "  ", -offset);
    vector<uint64_t> sigmas(N+1, 1);

    for (uint32_t i = start_index; i <= N; i++) {
        uint64_t num = (uint64_t) square_mult * i * i - offset;
        uint64_t rem = num / status[i].product;
        uint64_t sigma = status[i].sigma;

        // Handle 2's
        if (rem && (rem & 1) == 0) {
            uint32_t twos = 0;
            while (rem && (rem & 1) == 0) {
                rem >>= 1;
                twos += 1;
            }
            sigma *= (1 << (twos + 1)) - 1;
        }

        //printf("\t\t%u -> %lu | %lu & %lu\n", i, num, rem, sigma);

        // Handle any remaining prime
        if (rem > 1) {
            assert((rem > N) && (rem % 2 == 1));
            sigma *= (1 + rem);
        }
        sigmas[i] = sigma;
    }
    return sigmas;
}

int main(int argc, char** argv) {
    // Allow comma seperators
    setlocale(LC_NUMERIC, "");

    const uint64_t N = (argc != 2) ? 1e12 : atol(argv[1]);
    printf("Testing up to %'ld (%.1e)\n", N, (double) N);
    /**
     * First number with sigma > uint64_t
     * A002093(3491) = 3590449939146470400
     * A034885(3491) = sigma(A002093(3491)) = 18451979754511564800
     */
    if (N >= 3'590'449'939'146'470'400ul) {
        printf("sigma(i) can overflow with N=%.1e\n", (double) N);
        exit(1);
    }

    const uint32_t STOP = isqrt(N);
    const uint32_t HALF_STOP = isqrt(N/2);

    vector<uint32_t> primes = gen_primes(3, STOP);
    assert(primes.front() == 3);

    vector<uint32_t> half_primes;
    for (auto p : primes) { if (p <= HALF_STOP) half_primes.push_back(p); };

    printf("primes(%'u) = |%lu| %u %u %u ... %u\n",
            HALF_STOP, half_primes.size(), half_primes[0], half_primes[1], half_primes[2], half_primes.back());
    printf("primes(%'u) = |%lu| %u %u %u ... %u\n",
            STOP, primes.size(), primes[0], primes[1], primes[2], primes.back());
	printf("\n");

    /**
     * Compute both the square and non-square sequence first
     * This increases max memory 25% (2 + 2 + 1 vs 2+2)
     * BUT allows us to process both positive and negative offset at the same time
     */

    std::function factor_offset_bounded{
        [&](bool twice_square, int32_t offset) {
            if (twice_square) {
                return factor_offset(HALF_STOP, true, offset, half_primes);
            } else {
                return factor_offset(STOP, false, offset, primes);
            }
        }
    };

    const auto sigmas_squares = factor_offset_bounded(/* twice_square= */ false, /* offset= */ 0);
    const auto sigmas_twice_squares = factor_offset_bounded(/* twice_square= */ true, /* offset= */ 0);

    uint32_t max_index = 0;

    #pragma omp parallel for schedule(dynamic, 1)
    for (int32_t offset_abs = 1; offset_abs <= 99999; offset_abs += 2) {
        uint32_t matches = 0;

        for (bool twice_square : {false, true}) {
            const auto& compare_with = twice_square ? sigmas_twice_squares : sigmas_squares;

            for (auto offset : {offset_abs, -offset_abs}) {
                auto sigmas_offset = factor_offset_bounded(twice_square, offset);
                assert(compare_with.size() == sigmas_offset.size());

                {
                    uint32_t square_mult = (1 + twice_square);
                    const uint32_t start_index = offset >= 0 ? 0 : isqrt((2 - offset + square_mult - 1) / square_mult);
                    for (uint32_t i = start_index; i < compare_with.size(); i++) {
                        if (compare_with[i] - offset == sigmas_offset[i]) {
                            matches += 1;
                            if (offset != 7) {
                                if (i > max_index) {
                                    max_index = i;
                                    printf("NEW MAX INDEX: %u\n", i);
                                }
                            }

                            uint64_t num = (uint64_t) square_mult * i * i;
                            /*
                            printf("\nDIST=%u MATCH @ %s%u^2, %lu and %lu | %lu %-d vs %lu\n\n",
                                abs(offset), twice_square ? "2*" : "", i,
                                num, num - offset,
                                compare_with[i], -offset, sigmas_offset[i]);
                            // */
                            // Format used for README.md
                            printf("  * DIST=%u:\t(%lu, %lu), `%lu=%s%u^2`\n",
                                offset_abs,
                                std::min(num, num - offset), std::max(num, num - offset),
                                num, twice_square ? "2*" : "", i);
                        }
                    }
                }
            }
        }
        if (matches > 1 && !(offset_abs == 7 || offset_abs == 1097)) {
            printf("MULTIPLE MATCHES!\n");
            exit(0);
        }
    }
}
