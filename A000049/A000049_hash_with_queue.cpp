#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <queue>
#include <set>
#include <unordered_set>
#include <utility>
#include <vector>

using std::pair;
using std::vector;
using std::priority_queue;

using std::cout;
using std::endl;

/**
 * Population of 3 x^2 + 4 y^2
 */

// {x, y}
typedef vector<std::pair<uint32_t, uint32_t>> congruence;


struct data
{
    uint64_t n_3x2p4y2;
    uint32_t x;
    uint32_t y;
};

class compGT
{
    public:
        bool operator() (const data& A, const data& B) const
        {
            return A.n_3x2p4y2 > B.n_3x2p4y2;
        }
};


/**
 * Expand one congruence class of the population
 *
 * Each (x, y) pair should have same n % base
 */
pair<uint64_t, uint64_t>
expand_class(uint64_t N, uint32_t base, congruence &parts) {

    uint64_t three_base_squared = 3ul * base * base;
    uint64_t four_base_squared = 4ul * base * base;
    uint64_t six_base = 6ul * base;
    uint64_t eight_base = 8ul * base;

    uint64_t population = 0;
    uint64_t enumerated = 0;

    priority_queue<data, std::vector<data>, compGT> items;
    data item;

    for (const auto& d : parts) {
        item.x = d.first;
        item.y = d.second;
        item.n_3x2p4y2 = 3ul * item.x * item.x + 4ul * item.y * item.y;
        assert(item.n_3x2p4y2 <= N);
        items.push(item);
    }

    uint32_t is_first = (items.top().n_3x2p4y2 % base) == 1;
    if (is_first) {
        cout << "\tpriority_queue start size: " << items.size() << endl;
    }

    // Extra item so don't have to check for empty
    items.push({N+1, 0, 0});

    uint64_t last_n = N+1; // so (0,0) doesn't match
    while (items.top().n_3x2p4y2 <= N)
    {
        item = items.top();
        items.pop();

        if (item.y < base) {
            // add {x + base, y}
            uint64_t n = item.n_3x2p4y2 + 6ul * base * item.x + three_base_squared;
            items.push({n, item.x + base, item.y});
        }

        enumerated++;
        population += item.n_3x2p4y2 != last_n;

        if (item.n_3x2p4y2 == last_n)
        {
            // Increment all items with same n
            while (items.top().n_3x2p4y2 == last_n)
            {
                enumerated++;
                data tempItem = items.top();
                items.pop();

                if (tempItem.y < base) {
                    // add {x + base, y}
                    uint64_t n = tempItem.n_3x2p4y2 + six_base * tempItem.x + three_base_squared;
                    if (n <= N) {
                        items.push({n, tempItem.x + base, tempItem.y});
                    }
                }

                tempItem.n_3x2p4y2 += eight_base * tempItem.y + four_base_squared;
                tempItem.y += base;
                if (tempItem.n_3x2p4y2 <= N)
                    items.push(tempItem);
            }
        }

        last_n = item.n_3x2p4y2;

        // 4*((y + base)^2 - y^2) = 4 * (2*base*y + base*base)
        item.n_3x2p4y2 += eight_base * item.y + four_base_squared;
        item.y += base;
        if (item.n_3x2p4y2 <= N)
            items.push(item);
    }
    //cout << "\t" << population << " / " << enumerated << endl;

    return {enumerated, population};
}

vector<congruence> build_congruences(uint64_t N, uint64_t num_classes)
{
    // Quit if not enough memory (~8GB) to store all congruence classes.
    if ((num_classes * num_classes * 9) > (1ull << 33)) {
        fprintf(stderr, "TOO MANY CLASES %lu\n", num_classes);
        exit(1);
    }

    vector<congruence> classes(num_classes);
    classes[0].reserve(2 * num_classes);
    for (uint32_t r = 1; r < num_classes; r++) {
        classes[r].reserve(num_classes+1);
    }

    uint64_t elements = 0;
    for (uint32_t x = 0; x < num_classes ; x++) {
        uint64_t temp_x = (uint64_t) 3ul * x * x;
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
    const uint64_t num_classes = 1009;

    vector<congruence> classes = build_congruences(N, num_classes);

    uint64_t elements = 0;
    for (size_t cls = 0; cls < num_classes; cls++)
        elements += classes[cls].size();

    setlocale(LC_NUMERIC, "");
    printf("\tnum_clases: %lu\n", num_classes);
    printf("\telements: %'lu (<= %lu)\n", elements, num_classes * num_classes);
    printf("\tmemory: ~%.1f MB\n", 9.0 * elements / 1024 / 1024);
    assert(elements <= (num_classes * num_classes));
    float guess_pop_per = (float) N / (15 - bits / 5) / num_classes;
    printf("\tpopulation per residual ~%.0f\n", guess_pop_per);

    uint64_t population = 0;
    uint64_t enumerated = 0;

    const uint64_t CPU_SPLIT = 128;
    // Outer loop to parallel without contention on Set
    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t v = 0; v < CPU_SPLIT; v++) {
        for (size_t m = v; m < num_classes; m += CPU_SPLIT) {

            auto [enumerated_class, population_class] = expand_class(N, num_classes, classes[m]);

            #pragma omp critical
            {
                enumerated += enumerated_class;
                population += population_class;
            }
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
