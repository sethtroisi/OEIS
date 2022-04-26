#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <queue>
#include <set>
#include <unordered_set>
#include <utility>
#include <vector>

using std::vector;

using std::cout;
using std::endl;

/**
 * Population of 3 x^2 + 4 y^2
 */

// {x, y}
typedef vector<std::pair<uint32_t, uint32_t>> congruence;

//typedef std::unordered_set<uint64_t> Set;
//typedef std::unordered_set<uint32_t> Set;
 #include "flat_hash_map.hpp"
 typedef   ska::flat_hash_set<uint32_t> Set;

/**
 * Expand one congruence class of the population
 *
 * Each (x, y) pair should have n % base == residual
 */
std::pair<uint64_t, uint64_t>
expand_class(
        uint64_t N, uint64_t mod_base, uint64_t residual,
        size_t num_passes,
        Set &found,
        congruence &parts) {

    size_t shift = 0;
    while (2ul << (++shift) <= mod_base);
    // X >> shift is unique if X is a multiple of mod_base
    assert((1ul << shift) <= mod_base);
    assert((2ul << shift) > mod_base);

    uint64_t four_base_squared = 4ul * mod_base * mod_base;
    uint64_t eight_base = 8ul * mod_base;
    uint64_t eight_base_squared = 8ul * mod_base * mod_base;

    // Do several passes over the data
    // On each pass
    //   * add a few values for each x (up to a new min)
    //   * Do a removal pass over found (removing everything below new min)
    //      * Should be >75%

    // Build list of all (x, y)
    vector<std::pair<uint64_t, uint32_t>> X;
    {
        for (const auto& d : parts) {
            for (uint32_t x = d.first; ; x += mod_base) {
                uint64_t temp_x = 3ul * x * x;
                if (temp_x > N)
                    break;

                X.push_back({temp_x, d.second});
            }
        }
        // To slow to be valueable
        //std::sort(X.begin(), X.end());
        if (residual == 1)
            printf("\t%lu X values\n", X.size());

        parts.clear();
    }

    // Max item that can be inserted/added to found
    const size_t found_max_element = std::numeric_limits<
        typename std::remove_reference_t<decltype(found)>::key_type
    >::max();

    uint64_t found_count = 0;
    uint64_t enumerated = 0;

    for (size_t pass = 0; pass < num_passes; pass++) {
        // Count number of values [pass_min, pass_max];
        size_t pass_min = (__uint128_t) N * pass / num_passes + 1;
        size_t pass_max = (__uint128_t) N * (pass + 1) / num_passes;
        if (pass == 0) {
            pass_min = 0;
        }

        size_t pass_enumerated = 0;

        const size_t max_element = (pass_max - pass_min) >> shift;
        assert(max_element <= found_max_element);

        for(auto& d : X) {
            /* 3*x^2 */
            uint64_t temp_x = d.first;
  //          if (temp_x > pass_max) Requires sort(X) which is slow
  //              break;

            uint32_t y = d.second;

            // 4 * ((y + base)^2 - y^2) = 8*base*y + 4*base^2
            // derivative with y and y+base => 8*base*base
            uint64_t temp_y = ((uint64_t) y * y) << 2;

            uint64_t n = temp_x + temp_y;
            uint64_t y_delta = eight_base * y + four_base_squared;

            for (; n <= pass_max;) {
                //assert(n % mod_base == residual);
//                assert(n >= pass_min);
//                assert(n <= pass_max);

//                found.insert(n);
                found.insert((n - pass_min) >> shift);
                pass_enumerated++;

                // N takes value for new y, but not inserted into found yet
                n += y_delta;
                y_delta += eight_base_squared;
                y += mod_base;
            }

            // Save progress to y
            d.second = y;
        }

        if (residual == 1)
            printf("\tpass %2lu [%lu, %lu] -> %lu/%lu\n",
                    pass, pass_min, pass_max,
                    found.size(), pass_enumerated);

        found_count += found.size();
        enumerated += pass_enumerated;
        found.clear();
    }

    return {found_count, enumerated};
}

vector<congruence> build_congruences(uint64_t N, uint64_t num_classes)
{
    // Quit if not enough memory (~8GB) to store all congruence classes.
    if ((num_classes * num_classes * 9) > (1ull << 33)) {
        fprintf(stderr, "TOO MANY CLASSES %lu\n", num_classes);
        exit(1);
    }

    vector<congruence> classes(num_classes);
    classes[0].reserve(2 * num_classes);
    for (uint32_t r = 1; r < num_classes; r++) {
        classes[r].reserve(num_classes+1);
    }

    uint64_t elements = 0;
    for (uint32_t x = 0; x < num_classes ; x++) {
        uint64_t temp_x = (uint64_t) x * x;
        temp_x += temp_x << 1;
        if (temp_x > N)
            break;

        uint64_t temp_n = temp_x;
        // 4 * (y + 1) ^ 2 = 4 * y^2 + 8*y + 4;
        uint32_t delta_y = 4;

        for (uint32_t y = 0; y < num_classes && temp_n < N; y++) {
            elements++;

            uint32_t cls = temp_n % num_classes;
            classes[cls].emplace_back(x, y);

            temp_n += delta_y;
            delta_y += 8;
        }
    }

    for (size_t cls = 0; cls < num_classes; cls++) {
        //fprintf(stderr, "\t%lu -> %lu\n", cls, classes[cls].size());
        if (cls > 0)
            assert(classes[cls].size() <= num_classes + 1);
    }

    return classes;
}


int main(int argc, char** argv)
{
    auto start = std::chrono::steady_clock::now();

    size_t bits = 25;
    if (argc == 2) {
        bits = atoi(argv[1]);
    }

    uint64_t N = 1ull << bits;


    // All congruence classes are only possible if num_classes is a prime
    //      4*k + 1 -> quadratic residual -> twice as many entries for 0
    //      4*k + 3 -> none quad residual -> 1 entry for 0
    // 37, 101, 331, 1009, 3343, 10007, 30011
    uint64_t num_classes = 10007; //30011;


    vector<congruence> classes = build_congruences(N, num_classes);

    uint64_t elements = 0;
    for (size_t cls = 0; cls < num_classes; cls++)
        elements += classes[cls].size();

    setlocale(LC_NUMERIC, "");
    printf("\tnum_classes: %lu\n", num_classes);
    printf("\telements: %'lu (<= %lu)\n", elements, num_classes * num_classes);
    printf("\tmemory: ~%.1f MB\n", 9.0 * elements / 1024 / 1024);
    assert(elements <= (num_classes * num_classes));
    float guess_pop_per = (float) N / (14 - bits / 9) / num_classes;
    printf("\tpopulation per residual ~%.0f\n", guess_pop_per);


    /**
     * Large values reduce Hash size BUT increase number of iterations over X
     * Aim for a 0.5-2 values of y per pass?
     */
    float X_per = sqrt(N / 3.0);
    size_t num_passes = 1 + 2 * guess_pop_per / X_per;
    printf("\tX_per ~%.0f -> %lu passes\n", X_per, num_passes);

    uint64_t population = 0;
    uint64_t enumerated = 0;

    const uint64_t CPU_SPLIT = 128;
    // Outer loop to parallel without contention on Set
    //#pragma omp parallel for schedule(dynamic, 1)
    for (size_t v = 0; v < CPU_SPLIT; v++) {
        // Allocated once here to avoid lots of memory allocation.
        Set found;
        //found.reserve(guess_pop_per / 0.3);

        uint64_t iter_cpu = 0;
        auto start_cpu = std::chrono::steady_clock::now();

        for (size_t m = v; m < num_classes; m += CPU_SPLIT) {
            found.clear();

            auto [f_class, e_class] = expand_class(
                N, num_classes, m, num_passes,
                found,
                classes[m]);

            //#pragma omp critical
            {
                population += f_class;
                enumerated += e_class;
                iter_cpu += e_class;
            }
        }

        if (v <= 10) {
            auto end = std::chrono::steady_clock::now();
            double elapsed = std::chrono::duration<double>(end-start_cpu).count();
            printf("\t\t%2lu, iters: %-12lu  iter/s: %.2f million\n",
                v, iter_cpu, iter_cpu / 1e6 / elapsed);
        }

    }

    // 0 doesn't count for this sequence.
    population -= 1;

    auto end = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(end-start).count();
    printf("| %2lu | %-12lu | %-12lu | %.1f | unique: %.2f  iter/s: %.1f million\n",
        bits, population, enumerated,
        elapsed, (float) population / enumerated, enumerated / 1e6 / elapsed);

}
