#include <cstdint>
#include <iostream>

#include "../utils/count_special_primes.hpp"

using std::cout;
using std::endl;


int main(int argc, char** argv) {
    size_t bits = 30;
    if (argc == 2) {
        bits = atoi(argv[1]);
    }

    const auto count = count_population_quadratic_form(
        bits,
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

    cout << "A000047(" << bits << ") = " << count << endl;
    return 0;
}
