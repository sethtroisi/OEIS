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
        /* startprime= */ 5,
        /* add_to_special_primes= */ 2,
        [](uint64_t n) {
            uint64_t m = n % 6;
            return (n / 6) + (m >= 1) +   (n >= 3);
        },
        [](uint64_t n) {
            uint64_t m = n % 6;
            return (n / 6) + (m >= 5) +   (n >= 2);
        },
        [](uint64_t p) { return (p % 6 == 1); }
    );

    cout << "A000205(" << bits << ") = " << count << endl;
    return 0;
}
