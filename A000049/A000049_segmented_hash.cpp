#include <algorithm>
#include <bitset>
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

/* Per Core Cache
 * Xeon W-2135 Optimal ~ 4*1024*1024 which is 50% of L2, 33% of L3
 * Ryzen 3900x Optimal ~ 32*1024*1024 which is 50% of L3 (over 8 threads)
 *
 * Larger caches significantly lower num_passes but have 1/10 the write speed
 * Which never makes up for reduced overhead.
 */
typedef std::bitset<32 * 1024 * 1024 + 1> Set;


/**
 * Expand one congruence class of the population
 *
 * Each (x, y) pair should have n % base == residual
 */
std::pair<uint64_t, uint64_t>
expand_class(uint64_t N, uint64_t mod_base, uint64_t residual, congruence &parts) {
    Set found;

    size_t shift = 0;
    {
      // X >> shift is unique if X is a multiple of mod_base
      while (2ul << (++shift) <= mod_base);
      assert((1ul << shift) <= mod_base);
      assert((2ul << shift) > mod_base);
    }

    size_t num_passes = ((N >> shift) - 1) / found.size() + 1;
    if (residual == 1) {
        printf("\tbitset<%lu> -> %lu passes\n", found.size(), num_passes);
    }
    // Needed for which pass pair is first included in, should be slightly smaller than IRL
    uint64_t size_per_pass = N / num_passes + 1;

    uint64_t four_base_squared = (uint64_t) 4ul * mod_base * mod_base;
    uint64_t eight_base_squared = 2ul * four_base_squared;
    uint64_t eight_base = 8ul * mod_base;

    // Build list of all (3*x^2, y_delta)
    // y_delta can almost be uint32_t but breaks eventually
    vector<vector<std::pair<uint64_t, uint64_t>>> X;
    X.resize(num_passes);
    {
        for (const auto& [x1, y] : parts) {
            uint64_t temp_y = 4ul * y * y;
            // 4 * ((y + base)^2 - y^2) = 8*base*y + 4*base^2
            uint64_t y_delta = eight_base * y + four_base_squared;

            // (0,0) -> 0 isn't "valid";
            uint32_t x = (x1 == 0 && y == 0) ? mod_base : x1;
            for (; ; x += mod_base) {
                uint64_t temp_n = 3ul * x * x + temp_y;
                if (temp_n > N)
                    break;

                // Pseudo radix sort! Determines the first pass that needs (x, y)
                // This can underestimate by one to ease math requirement
                uint32_t first_pass = temp_n / size_per_pass;

                assert( temp_n >= (__uint128_t) N * first_pass / num_passes + 1 );
                assert(0 <= first_pass);
                assert(first_pass < num_passes);
                X[first_pass].push_back({temp_n, y_delta});
            }
        }
        parts.clear();
        if (residual == 1) {
            size_t num_X = 0;
            for (const auto& t : X) num_X += t.size();
            printf("\tresidual %ld |pairs| = %lu/%lu\n", residual, num_X, num_passes);
        }
    }

    auto start_class = std::chrono::high_resolution_clock::now();
    uint64_t total_found = 0;
    uint64_t total_enumerated = 0;

    for (size_t pass = 0; pass < num_passes; pass++) {
        // Count number of values [pass_min, pass_max];
        size_t pass_min = (__uint128_t) N * pass / num_passes + 1;
        size_t pass_max = (__uint128_t) N * (pass + 1) / num_passes;

        // Numbers included in interval (+1 as both endpoints are included)
        //size_t pass_interval_size = pass_max - pass_min + 1;
        const size_t max_element = (pass_max - pass_min) >> shift;
        assert(max_element < found.size());

        size_t pass_enumerated = 0;
        size_t pass_iterated = 0;
        for(size_t i = 0; i <= pass; i++) {
            pass_iterated += X[i].size();
            for(auto& d : X[i]) {
                uint64_t n = d.first;
                uint64_t y_delta = d.second;

                for (; n <= pass_max;) {
                //for (; n <= pass_interval_size;)
                    //assert(n % mod_base == residual);
                    assert(n >= pass_min);
                    assert(n <= pass_max);
                    //assert(n <= pass_interval_size);

                    found.set((n - pass_min) >> shift);
                    //found.set(n >> shift);
                    pass_enumerated++;

                    // N takes value for new y, but not inserted into found yet
                    n += y_delta;
                    y_delta += eight_base_squared;
                }

                // Subtract off interval, to avoid subtraction in the inner loop
                // assert (n > pass_interval_size);
                // n -= pass_interval_size;

                // Save ending point of this pass (starting point of next pass)
                d.first = n;
                d.second = y_delta;
            }
        }

        size_t pass_found = found.count();

        if (residual == 1 && (
                    (pass + 1 == num_passes) ||
                    (pass <= 4) ||
                    (pass <= 128 && pass % 16 == 0) ||
                    (pass % 128 == 0))) {
            printf("\tpass %2lu [%lu, %lu] -> %lu/%lu/%lu\n",
                    pass, pass_min, pass_max,
                    pass_found, pass_enumerated, pass_iterated);
        }

        total_found += pass_found;
        total_enumerated += pass_enumerated;

        if (pass != num_passes - 1) {
          found.reset();
        }
    }

    if (residual == 1 || N >= (1ul << 46)) {
        auto end = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double>(end - start_class).count();
        printf("\tresidual %lu, iters: %lu secs: %.2f -> %.1f million iter/s\n",
          residual, total_enumerated, elapsed, total_enumerated / 1e6 / elapsed);
    }
    return {total_found, total_enumerated};
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

        for (uint32_t y = 0; y < num_classes && temp_n <= N; y++) {
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
    auto start = std::chrono::high_resolution_clock::now();

    size_t bits = 25;
    if (argc == 2) {
        bits = atoi(argv[1]);
    }

    uint64_t N = 1ull << bits;

    /*
     *  All congruence classes are only possible if num_classes is a prime
     *      4*k + 1 -> quadratic residual -> twice as many entries for 0
     *      4*k + 3 -> none quad residual -> 1 entry for 0
     *
     * For the bitset approach it's best to set the smallest number that doesn't
     * explode num_passes
     */
    // 37, 101, 331, 1031, 2053, 4099, 8209, 16411, 32771
    uint64_t num_classes = 2053; //4099; //8209;

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

    uint64_t population = 0;
    uint64_t enumerated = 0;

    const uint64_t CPU_SPLIT = 128;
    // Outer loop to parallel without contention on Set
    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t v = 0; v < CPU_SPLIT; v++) {
        // v = 0 is weird, do it last
        size_t w = (v + 1) % CPU_SPLIT;

        uint64_t iter_cpu = 0;
        auto start_cpu = std::chrono::high_resolution_clock::now();

        for (size_t c = w; c < num_classes; c += CPU_SPLIT) {
            auto [f_class, e_class] = expand_class(N, num_classes, c, classes[c]);

            #pragma omp critical
            {
                population += f_class;
                enumerated += e_class;
                iter_cpu += e_class;
            }
        }

        // I wish #pragma ordered wasn't broken
        if (w <= 4 || w == 8 || w == 16 || (w % 32 == 0)) {
            auto end = std::chrono::high_resolution_clock::now();
            double elapsed = std::chrono::duration<double>(end-start_cpu).count();
            double total_elapsed = std::chrono::duration<double>(end-start).count();
            printf("\t\t%2lu, iters: %-12lu  iter/s: %.2f / %.1f\n",
                w,
                iter_cpu,
                iter_cpu / 1e6 / elapsed,
                enumerated / 1e6 / total_elapsed);
        }

    }

    auto end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(end-start).count();
    printf("| %2lu | %-12lu | %-12lu | %.1f | iter/s: %.1f million\n",
        bits, population, enumerated, elapsed, enumerated / 1e6 / elapsed);
}
