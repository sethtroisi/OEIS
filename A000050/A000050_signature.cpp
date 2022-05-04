#include <algorithm>
#include <cassert>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <vector>

#include <primesieve.hpp>
#include "../utils/count_special_primes.hpp"


using std::vector;

using std::cout;
using std::cerr;
using std::endl;


int main(int argc, char** argv) {
    size_t bits = 30;
    if (argc == 2) {
        bits = atoi(argv[1]);
    }

    const auto count = count_population_quadratic_form(
        bits,
        /* start_prime= */ 3,
        /* add_to_special_primes= */ 0,
        [](uint64_t n) { return (n / 4) + ((n % 4) >= 1); },
        [](uint64_t n) { return (n / 4) + ((n % 4) >= 3); },
        [](uint64_t p) { return (p & 3) == 1; }
    );

    cout << "A000050(" << bits << ") = " << count << endl;
    return 0;
}
