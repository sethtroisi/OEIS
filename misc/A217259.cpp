// g++ --std=c++11 -O3 -Werror -Wall A217259.cpp -lgmpxx -lgmp

#include <cassert>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <vector>

#include <gmp.h>
#include <gmpxx.h>

using std::vector;

/**
 * Calculate sigma[i] for i in range(start, start+n).
 */
vector<uint64_t> SegmentedSieveOfSigma(uint64_t start, uint64_t N) {
    auto past = start + N;
    auto sigmas = vector<uint64_t>(N, 1);

    // HACK for sigmas(0) and sigmas(1)
    for (int i = 0; i < 2; i++) {
        if (start + i < 2) {
            sigmas[i] = 0;
        }
    }

    // Adjust to include n as a divisor of n
    for (uint64_t i = 0; i < N; i++) {
        sigmas[i] += start + i;
    }

    // isqrt with gmp
    mpz_class sqrt = past-1;
    mpz_sqrt(sqrt.get_mpz_t(), sqrt.get_mpz_t());
    uint64_t isqrt = mpz_get_ui(sqrt.get_mpz_t());
    assert( isqrt * isqrt < past );
    assert( (isqrt+1) * (isqrt+1) >= past );

    for (uint64_t factor = 2; factor <= isqrt; factor++) {
        auto f2 = factor * factor;
        assert(f2 < past);

        uint64_t ceil, next_index;

        if (f2 >= start) {
            assert(f2 - start < N);
            // n=factor^2 only gets to count factor, not factor + n/factor
            sigmas[f2 - start] += factor;

            ceil = factor + 1;
            next_index = f2 + factor - start;
        } else {
            // ceil(start / factor)
            ceil = (start-1) / factor + 1;
            next_index = ceil * factor - start;
        }

        uint64_t count = ceil;
        for (uint64_t index = next_index; index < N; index += factor, count++) {
            // count = number / factor
            sigmas[index] += factor + count;
        }
    }

    return sigmas;
}



//S = time.time()
void print_match(uint64_t mid) {
    static uint32_t found = 0;
    static uint64_t print_mult = 1;
    static auto S = std::chrono::system_clock::now();
    static auto next_time = 5;

    std::chrono::duration<double> elapsed = std::chrono::system_clock::now() - S;

    // TODO commas in mid

    found += 1;
    if (found % print_mult == 0) {
        printf("%-8d %-'13lu\n", found, mid);
        if (found == 10 * print_mult)
            print_mult *= 10;
    } else if (elapsed.count() > next_time) {
        float rate = mid / elapsed.count() / 1e6;
        printf("%-8d %-'13lu\t\t%.1f seconds elapsed %.1fM/s\n",
                found, mid, elapsed.count(), rate);
        next_time += 5;
    }

    if (mid == 435 or mid == 8576 or mid == 8826)
        return;

    // Verify mid-1 and mid+1 are prime
    mpz_class mid_m_1 = mid - 1;
    mpz_class mid_p_1 = mid + 1;
    if((mpz_probab_prime_p(mid_m_1.get_mpz_t(), 25) != 2) ||
       (mpz_probab_prime_p(mid_p_1.get_mpz_t(), 25) != 2)) {
        printf("WHAT %lu | %lu\n", mid - 1, mid + 1);
        exit(1);
    }
}

void iterate(uint64_t START, uint64_t STOP, uint64_t SEGMENT) {
    // sigma(start-2), sigma(start-1)
    uint64_t last_sigmas[2] = {0, 0};

    if (START > 0)
        printf("\tCAN'T CHECK %lu or %lu\n", START, START+1);

    for (uint64_t start = START; start <= STOP; start += SEGMENT) {
        auto sigmas = SegmentedSieveOfSigma(start, SEGMENT);

        if (sigmas[0] - last_sigmas[0] == 2)
            print_match(start - 1);

        if (sigmas[1] - last_sigmas[1] == 2)
            print_match(start);

        for (uint32_t i = 1; i < SEGMENT-1; i++) {
            if (sigmas[i+1] - sigmas[i-1] == 2) {
                print_match(start + i);
            }
        }

        last_sigmas[0] = sigmas[SEGMENT-2];
        last_sigmas[1] = sigmas[SEGMENT-1];
    }
}


int main() {
    // Allow comma seperators
    setlocale(LC_NUMERIC, "");

    uint64_t START = 0;
    uint64_t SEGMENT = 1 << 17;
    uint64_t STOP = 1e15;
    iterate(START, STOP, SEGMENT);
}

// N < 10^8 | 5.3 seconds w/ SEGMENT = 2^16 | 440315 99999588
// N < 10^9 | 48  seconds w/ SEGMENT = 2^16 | 3424509 999999192
// N < 10^9 | 56  seconds w/ SEGMENT = 2^20 | 3424509 999999192
