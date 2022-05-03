#include <cassert>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <vector>

#include <primecount.hpp>
#include <primesieve.hpp>
#include "count_special_primes.hpp"


using std::vector;

using std::cout;
using std::cerr;
using std::endl;


void test_primepi_speed(size_t bits) {
    uint64_t n = 1ul << bits;
    uint64_t r = sqrt(n) + 1;
    while (r*r > n) {
        r--;
    }
    assert(r*r <= n);
    assert((r+1) * (r+1) > n);
    assert(r < std::numeric_limits<uint32_t>::max());

    auto start = std::chrono::high_resolution_clock::now();

    const auto count_special_primes = get_special_prime_counts(
        n, r,
        /* start_prime= */ 3,
        [](uint64_t n) { return (n / 4) + ((n % 4) >= 1); },
        [](uint64_t n) { return (n / 4) + ((n % 4) >= 3); },
        [](uint64_t p) { return (p & 3) == 1; }
    );

    auto end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(end - start).count();
    fprintf(stderr, "\touter count_special_primes(%lu) = %lu  (%.1f)\n",
        n, count_special_primes.at(n), elapsed);

    primecount::set_num_threads(4);

    start = std::chrono::high_resolution_clock::now();

    vector<uint64_t> k;
    for (int i = 1; i < r; i++) {
        uint64_t v = n / i;
        k.push_back(primecount::pi(v));
    }

    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration<double>(end - start).count();
    fprintf(stderr, "\tcount_special_primes(%lu) = %lu  (%.1f)\n",
        n, count_special_primes.at(n), elapsed);
}

int main(int argc, char** argv) {
    size_t bits = 30;
    if (argc == 2) {
        bits = atoi(argv[1]);
    }

    test_primepi_speed(bits);
    return 0;
}
