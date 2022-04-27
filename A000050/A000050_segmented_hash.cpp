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

/**
 * Population of x^2 + y^2
 * (x <= y)
 */

// {x, y}
typedef vector<std::pair<uint32_t, uint32_t>> congruence;

/* Per Core Cache
 * Xeon W-2135 Optimal ~ 4*1024*1024 which is 50% of L2, 33% of L3
 * Ryzen 3900x Optimal ~ 32*1024*1024 which is 50% of L3 (over 8 threads)
 */
typedef   std::bitset<32 * 1024 * 1024 + 1> Set;


/**
 * Expand one congruence class of the population
 *
 * Each (x, y) pair should have n % base == residual
 */
std::pair<uint64_t, uint64_t>
expand_class(
        uint64_t N, uint64_t mod_base, uint64_t residual,
        congruence &parts) {
    Set found;

    size_t shift = 0;
    {
      // X >> shift is unique if X is a multiple of mod_base
      while (2ul << (++shift) <= mod_base);
      assert((1ul << shift) <= mod_base);
      assert((2ul << shift) > mod_base);
    }
    assert((N >> shift) > 0);

    size_t num_passes = ((N >> shift) - 1) / found.size() + 1;
    if (residual == 0) {
        printf("\tbitset<%lu> -> %lu passes\n", found.size(), num_passes);
    }

    uint64_t base_squared = (uint64_t) mod_base * mod_base;
    uint64_t two_base_squared = 2ul * base_squared;
    uint64_t four_base_squared = 4ul * base_squared;
    uint64_t two_base = 2ul * mod_base;

    // Build list of all (x^2 + y^2, y_delta)
    vector<std::pair<uint64_t, uint64_t>> X;
    {
        for (const auto& [x1, y1] : parts) {
            uint64_t y = y1;
            uint64_t x = x1;

            if (y < x) {
                y += mod_base;
            }
            assert(x <= y);
            assert(x + mod_base > y);

            uint64_t x_2 = x * x;
            uint64_t y_2 = y * y;
            // ((y + base)^2 - y^2) = 2*base*y + base^2
            uint64_t y_delta = two_base * y + base_squared;
            uint64_t n_delta = two_base * (x + y) + two_base_squared;

            uint64_t temp_n = x_2 + y_2;
            for (; temp_n <= N; ) {
                // TODO instead of pushing to X push
                // push to X_i where i is the first pass this will be included in
                //assert(x <= y);
                X.push_back({temp_n, y_delta});

                temp_n += n_delta;
                n_delta += four_base_squared;
                y_delta += two_base_squared;
                //x += mod_base;
                //y += mod_base;  // Needed to maintain the invariant x <= y
                //assert(temp_n == (x*x + y*y));
            }
        }
        if (residual <= 1)
            printf("\tresidual %ld |pairs| = %lu\n", residual, X.size());

        parts.clear();
    }

    auto start_class = std::chrono::high_resolution_clock::now();
    uint64_t total_found = 0;
    uint64_t total_enumerated = 0;

    for (size_t pass = 0; pass < num_passes; pass++) {
        // Count number of values [pass_min, pass_max];
        size_t pass_min = (__uint128_t) N * pass / num_passes + 1;
        size_t pass_max = (__uint128_t) N * (pass + 1) / num_passes;

        const size_t max_element = (pass_max - pass_min) >> shift;
        assert(max_element < found.size());

        // Hack for (0,0)
        if (residual == 0 && pass == 0) {
            auto& t = X[0];
            assert(t.first == 0);
            t.first += t.second;
            t.second += two_base_squared;
        }

        size_t pass_enumerated = 0;
        for(auto& d : X) {
            uint64_t n = d.first;
            uint64_t y_delta = d.second;

            for (; n <= pass_max;) {
                //assert(n == (x*x + y*y));
                //assert(n % mod_base == residual);
                //assert(n >= pass_min);
                //assert(n <= pass_max);

                found.set((n - pass_min) >> shift);
                //std::cout << n << std::endl;
                pass_enumerated++;

                // N takes value for new y, but not inserted into found yet
                n += y_delta;
                y_delta += two_base_squared;
                //y += mod_base;
            }

            // Save ending point of this pass (starting point of next pass)
            d.first = n;
            d.second = y_delta;
        }

        size_t pass_found = found.count();

        if (residual == 1 && (pass <= 4 || pass % 16 == 0)) {
            printf("\tpass %2lu [%lu, %lu] -> %lu/%lu\n",
                    pass, pass_min, pass_max,
                    pass_found, pass_enumerated);
        }

        total_found += pass_found;
        total_enumerated += pass_enumerated;

        if (pass != num_passes - 1) {
          found.reset();
        }
    }

    if (residual == 1) {
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

    // TODO only build with x <= y
    uint64_t elements = 0;
    for (uint64_t x = 0; x < num_classes ; x++) {
        uint64_t temp_x = x * x;
        if (temp_x > N)
            break;

        uint64_t temp_n = 2 * temp_x;
        // (y + 1)^2 - y^2 = 2*y + 1;
        uint32_t delta_y = 2 * x + 1;

        for (uint32_t y = x; y < num_classes && temp_n <= N; y++) {
            assert(temp_n == (x*x + y*y));
            elements++;

            uint32_t cls = temp_n % num_classes;
            classes[cls].emplace_back(x, y);
            // Note: a little silly to add to a vector verus duplicate item at far end
            if (x != y)
                classes[cls].emplace_back(y, x);

            temp_n += delta_y;
            delta_y += 2;
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
    uint64_t num_classes = 1031; //4099; //8209;

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
    //#pragma omp parallel for schedule(dynamic, 1)
    for (size_t v = 0; v < CPU_SPLIT; v++) {
        uint64_t iter_cpu = 0;
        auto start_cpu = std::chrono::high_resolution_clock::now();

        for (size_t m = v; m < num_classes; m += CPU_SPLIT) {
            auto [f_class, e_class] = expand_class(N, num_classes, m, classes[m]);

            //#pragma omp critical
            {
                population += f_class;
                enumerated += e_class;
                iter_cpu += e_class;
            }
        }

        // I wish #pragma ordered wasn't broken
        if (v <= 4 || v == 8 || v == 16 || (v % 32 == 0)) {
            auto end = std::chrono::high_resolution_clock::now();
            double elapsed = std::chrono::duration<double>(end-start_cpu).count();
            printf("\t\t%2lu, iters: %-12lu  iter/s: %.2f million\n",
                v, iter_cpu, iter_cpu / 1e6 / elapsed);
        }

    }

    auto end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(end-start).count();
    printf("| %2lu | %-13lu | %-13lu | %.1f | unique: %.2f  iter/s: %.1f million\n",
        bits, population, enumerated,
        elapsed, (float) population / enumerated, enumerated / 1e6 / elapsed);
}
