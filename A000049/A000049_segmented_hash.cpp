#include <bitset>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <limits>
#include <unordered_map>
#include <utility>
#include <vector>

using std::pair;
using std::unordered_map;
using std::vector;

using std::cout;
using std::endl;

/**
 * Population of 3 x^2 + 4 y^2
 */

// {x, y}
typedef vector<pair<uint32_t, uint32_t>> congruence;

/* Per Core Cache
 * Xeon W-2135 Optimal ~ 4*1024*1024 which is 50% of L2, 33% of L3
 * Ryzen 3900x Optimal ~ 32*1024*1024 which is 50% of L3 (over 8 threads)
 *
 * Larger caches significantly lower num_passes but have 1/10 the write speed
 * Which never makes up for reduced overhead.
 */
typedef std::bitset<32 * 1024 * 1024> Set;


/**
 * Expand one congruence class of the population
 *
 * Each (x, y) pair should have n % base == residual
 */
pair<uint64_t, uint64_t>
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

    // Verify adding 1 bit doesn't change num_passes needed
    assert( (((N >> shift) - 1) / (found.size()+1) + 1) == num_passes);

    // Needed for which pass pair is first included in, should be slightly smaller than IRL
    uint64_t size_per_pass = N / num_passes + 1;

    // For x_delta, y_delta respectively
    uint64_t three_base_squared = 3ul * mod_base * mod_base;
    uint64_t six_base_squared   = 6ul * mod_base * mod_base;
    uint64_t six_base           = 6ul * mod_base;
    uint64_t four_base_squared  = 4ul * mod_base * mod_base;
    uint64_t eight_base_squared = 8ul * mod_base * mod_base;
    uint64_t eight_base         = 8ul * mod_base;

    // Build list of all (3*x^2, y_delta)
    // y_delta can almost be uint32_t but breaks around 2^38
    vector<vector<pair<uint64_t, uint64_t>>> X;
    X.resize(num_passes);
    {
        for (const auto& [x_orig, y] : parts) {
            uint64_t x = x_orig;

            /**
             * (0,0) -> 0 isn't "a positive integer".
             * Add it's next child to the list (see logic in loop below)
             * And manually advance to the next item here.
             *
             * This is some Hacky code but I can't find a better solution.
             * this loop only runs a few million times so the overhead of
             * a 99.999% false condition here is much better than handling
             * this more natively in the enumeration loop below.
             */
            if (x_orig == 0 && y == 0) {
                uint64_t y = mod_base;
                uint64_t temp_y = 4ul * y * y;
                uint64_t y_delta = eight_base * y + four_base_squared;
                X[0].push_back({temp_y, y_delta});

                x += mod_base;
            }

            uint64_t temp_n = 3ul * x * x + 4ul * y * y;
            // (z + base)^2 - z^2 = 2*base * z + base^2
            uint64_t x_delta = six_base * x + three_base_squared;
            const uint64_t y_delta = eight_base * y + four_base_squared;

            assert( temp_n % mod_base == residual);
            // Trivially true but worth verifying
            assert( x_delta % mod_base == 0);
            assert( y_delta % mod_base == 0);

            for (; temp_n <= N; ) {
                // Pseudo radix sort! Determines the first pass that needs (x, y)
                // This can underestimate by one to ease math requirement
                uint32_t first_pass = temp_n / size_per_pass;
                assert( temp_n >= ((__uint128_t) N * first_pass / num_passes + 1) );
                assert( 0 <= first_pass && first_pass < num_passes);

                //assert( temp_n == 3ul * x * x + 4ul * y * y );
                //assert( temp_n % mod_base == residual);
                X[first_pass].push_back({temp_n, y_delta});

                temp_n += x_delta;

                x_delta += six_base_squared;
                x += mod_base;
            }

            // Verify something went wrong during incrementing
            assert( temp_n == (3ul * x * x + 4ul * y * y) );
            assert( x_delta == (six_base * x + three_base_squared) );
        }

        parts.clear();

        if (residual == 1) {
            size_t num_X = 0;
            for (const auto& t : X) num_X += t.size();
            printf("\tclass %-4ld |pairs| = %lu/%lu\n", residual, num_X, num_passes);
        }
    }

    uint64_t total_found = 0;
    uint64_t total_enumerated = 0;

    for (size_t pass = 0; pass < num_passes; pass++) {
        if (pass > 0) {
            found.reset();
        }

        // Count number of values [pass_min, pass_max];
        size_t pass_min = (__uint128_t) N * pass / num_passes + 1;
        size_t pass_max = (__uint128_t) N * (pass + 1) / num_passes;
        assert(pass_max <= N);

        const size_t max_element = (pass_max - pass_min) >> shift;
        assert(max_element < found.size());

        size_t pass_enumerated = 0;
        size_t pass_iterated = 0;
        for (size_t pass_i = 0; pass_i <= pass; pass_i++) {
            pass_iterated += X[pass_i].size();
            for(auto& d : X[pass_i]) {
                uint64_t n = d.first;
                uint64_t y_delta = d.second;
                assert(n >= pass_min);

                for (; n <= pass_max;) {
                    //assert(n % mod_base == residual);

                    found.set((n - pass_min) >> shift);
                    pass_enumerated++;

                    // N takes value for new y, but not inserted into found yet
                    n += y_delta;
                    y_delta += eight_base_squared;
                }

                // No need to assert(n > pass_max), definition of for-loop exit

                // Save ending point of this pass (starting point of next pass)
                d.first = n;
                d.second = y_delta;
            }
        }

        size_t pass_found = found.count();
        total_found += pass_found;
        total_enumerated += pass_enumerated;

        if (residual == 1 && (
                    (pass + 1 == num_passes) ||
                    (pass <= 4) ||
                    (pass <= 128 && pass % 16 == 0) ||
                    (pass % 128 == 0))) {
            printf("\t pass %4lu [%lu, %lu] -> %lu/%lu/%lu\n",
                    pass, pass_min, pass_max,
                    pass_found, pass_enumerated, pass_iterated);
        }
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
    assert((4ul * num_classes * num_classes) < std::numeric_limits<uint32_t>::max());

    vector<congruence> classes(num_classes);
    classes[0].reserve(2 * num_classes);
    for (uint32_t r = 1; r < num_classes; r++) {
        classes[r].reserve(num_classes+1);
    }

    for (uint32_t x = 0; x < num_classes ; x++) {
        uint64_t three_x_2 = 3ul * x * x;
        if (three_x_2 > N)
            break;

        // 4 * (y + 1) ^ 2 = 4 * y^2 + 8*y + 4;
        for (uint32_t y = 0; y < num_classes; y++) {
            uint64_t n = three_x_2 + 4ul * y * y;
            if (n > N)
                break;

            // TODO can reduce size of classes with
            // (x, y) -> (x, -y)

            uint32_t cls = n % num_classes;
            assert(0 <= cls && cls < num_classes);
            classes[cls].emplace_back(x, y);
        }
    }

    for (size_t cls = 0; cls < num_classes; cls++) {
        //fprintf(stderr, "\t%lu -> %lu\n", cls, classes[cls].size());
        if (cls > 0)
            assert(classes[cls].size() <= num_classes + 1);
    }

    return classes;
}


/**
 * Load finished classes from a temp file.
 * This is usefull for double checking, spot checking, and cloud instances
 */
unordered_map<uint32_t, pair<uint64_t, uint64_t>>
load_finished(std::string filename, size_t bits, uint64_t num_classes)
{
    unordered_map<uint32_t, pair<uint64_t, uint64_t>> loaded;
    int64_t N = 1ull << bits;

    std::ifstream file(filename);
    while (file.good()) {
        int64_t c, f, e;
        file >> c >> f >> e;
        if (!file.good())
            break;
        assert(0 <= c && (uint64_t) c < num_classes);
        assert(0 <= f && f < N);
        assert(1 <= e && e < N);

        assert(loaded.count(c) == 0);
        loaded[c] = {f, e};
    }
    file.close();

    if (loaded.size())
      cout << "\tLoaded " << loaded.size() << " pairs from \"" << filename << "\"" << endl;

    return loaded;
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
    uint64_t num_classes = 2053;

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

    std::string fn = "partial_" + std::to_string(bits)
        + "_" + std::to_string(num_classes) + ".txt";
    const auto finished_classes = load_finished(fn, bits, num_classes);
    std::ofstream outfile(fn, std::ios_base::app); // append

    // Helps stabalize iter/s. Otherwise build_congruences is slowly amortized.
    auto start2 = std::chrono::high_resolution_clock::now();

    uint64_t population = 0;
    uint64_t enumerated = 0;
    uint64_t run_enumerated = 0;

    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t c = 0; c < num_classes; c++) {
        auto start_class = std::chrono::high_resolution_clock::now();

        // 0 is weird so swap order of 0 and 1
        size_t class_i = c > 1 ? c : 1 - c;
        if (finished_classes.count(class_i)) {
            #pragma omp critical
            {
                const auto prev = finished_classes.at(class_i);
                population += prev.first;
                enumerated += prev.second;
            }
            continue;
        }

        const auto& [f_class, e_class] = expand_class(
            N, num_classes, class_i, classes[class_i]);

        #pragma omp critical
        {
            population += f_class;
            enumerated += e_class;
            run_enumerated += e_class;

            outfile << class_i << " " << f_class << " " << e_class << endl;
            //printf("%4lu -> %lu/%lu\n", class_i, f_class, e_class);

            if (c <= 4 || c == 8 || c == 16 ||
                    (c < 128 && c % 32 == 0) ||
                    (c % 128 == 0)) {
                auto end = std::chrono::high_resolution_clock::now();
                double elapsed = std::chrono::duration<double>(end - start_class).count();
                double total_elapsed = std::chrono::duration<double>(end - start2).count();
                printf("\tclass %-4lu, iters: %-12lu  iter/s: %.2f / %.1f\n",
                    class_i,
                    e_class,
                    e_class / 1e6 / elapsed,
                    run_enumerated / 1e6 / total_elapsed);
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(end-start).count();
    printf("| %2lu | %-12lu | %-12lu | %.1f | iter/s: %.1f million\n",
        bits, population, enumerated, elapsed, enumerated / 1e6 / elapsed);
}
