#include "count_special_primes.hpp"
#include <sys/types.h>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <functional>
#include <utility>
#include <vector>

#include <primesieve.hpp>

#define LIBDIVIDE_AVX2
#include <libdivide.h>

using std::pair;
using std::vector;

#define INNER_ASSERT 0

/**
 * TODO: IDEA from SQL join.
 *
 * In theory we do
 *      V[N/1] -= V[N/1/p_1] - V[p_1-1],
 *      V[N/2] -= V[N/2/p_1] - V[p_1-1],
 *      V[N/3] -= V[N/3/p_1] - V[p_1-1],
 *      ...
 *      V[p^2] -= V[p^2/p_1] - V[p_1-1],
 *
 * Then
 *
 *      V[N/1] -= V[N/1/p_1] - V[p_2-1],
 *      V[N/2] -= V[N/2/p_1] - V[p_2-1],
 *      V[N/3] -= V[N/3/p_1] - V[p_2-1],
 *      ...
 *      V[p^2] -= V[p^2/p_1] - V[p_2-1],
 *
 * What's safe to reorder?
 *    numbers larger than V[N/p_1] aren't used recursively
 *    numbers smaller than V[p_1^2] don't change
 *
 *    so the updates could be grouped
 *        V[N/i] -= (
 *            V[N/i/p_1] - V[p_1-1]
 *            V[N/i/p_2] - V[p_2-1]
 *        )
 *     This can maybe be rewritten as
 *        V[N/i] += V[p_1-1] + V[p_2-1] + V[p_3-1] + ...
 *                  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 *                  constant and doesn't depend on i!
 *        V[N/i] -= V[N/i/p_1] + V[N/i/p_2] + V[N/i/p_3] + ...
 *                  depends on i and requires a lot of gather
 *
 * 99.9%+ of all updates are to small numbers (e.g. V[k] where k <= r)
 *      These are harder to reorder because large V[N/i] can depend on them.
 * 99% of all updates are with primes < n^(1/3) so updates more updates affect other updates
 * 60% of all updates are with primes < n^(1/4)
 *
 * let V_p = the values of V after prime p
 *      V_p[:p^2] doesn't change
 *      V_pi depends on V[pi^2: N/pi]
 *         this is the "active" portion of V which is changing between each iteration
 *         need most intermediate values so all values have to be computed
 *         might still be some memory bandwidth advantage to chunking
 *
 * New plan to do less memory reads
 * for small primes (less than n^(1/4))
 *      lots of V[N/j/p] can be constant for many j
 *          Think about N = 10**6, r = 1000
 *              V has all values [1...1000] then 500 more values in [1000, 2000]
 *              For p = 31
 *                  V[496..527] are all changed by the same amount V[17] - V[32]
 *                      Thinking of something like bulk updating to V[p*i: p*i+p-1]
 *                          This kinda looks like the second half of the sqrt update trick.
 *
 * upper loop is N/j/p > r
 * lower loop is
 *      avx512 till floor(N/j/p) == floor(N/(j+10)/p)
 *          N/j/p = N/(j+10)/p + 1
 *          N/j = N/(j+10) + p
 *          N*(j+10) - N*j = p * j * (j+10)
 *          N*10 = p * j * (j+10)
 *          j ~= sqrt(N * 10 / p)
 *
 * IRL
 *      N = 2^48
 *      r = 2^24
 *      P = 137
 *          upper loop is j < 122,500
 *          middle is j < 4,532,727
 *          chunked is j < 16,777,216
 *          then chunks of exactly p size for 16,777,216 - 18769 items
 */

/** Setup vectors */
void inline setup_vectors(
    const uint64_t n, const uint32_t r,
    std::function< uint64_t(uint64_t)> init_count_a,
    std::function< uint64_t(uint64_t)> init_count_b,
    vector<uint64_t> &counts_backing_v,
    vector<uint64_t> &counts_backing_a,
    vector<uint64_t> &counts_backing_b
)
{
    const auto length = ((n/r)-1) + r;
    {
        counts_backing_v.reserve(length);
        counts_backing_a.reserve(length);
        counts_backing_b.reserve(length);
        // 1, 2, ... n / r - 1    <- r is always guaranteed to be in this sequence
        for(uint32_t v = 1; v < (n / r); v++) {
            uint64_t c_a = init_count_a(v);
            uint64_t c_b = init_count_b(v);
            assert(c_a + c_b <= v);
            counts_backing_v.push_back(v);
            counts_backing_a.push_back(c_a);
            counts_backing_b.push_back(c_b);
        }

        // n/r, n/(r-1), n/(r-2), ... n/3, n/2, n/1
        for(uint64_t i = r; i >= 1; i--) {
            uint64_t v = n / i;
            uint64_t c_a = init_count_a(v);
            uint64_t c_b = init_count_b(v);
            assert(c_a + c_b <= v);
            counts_backing_v.push_back(v);
            counts_backing_a.push_back(c_a);
            counts_backing_b.push_back(c_b);
        }
    }
}

/**
 * Process i = [i_bot, i_top] in reverse order
 */
__attribute__((always_inline)) inline
void middle_loop_avx(
    uint32_t i_bot, uint32_t i_top,
    uint64_t c_a, uint64_t c_b,
    uint64_t (&counts)[6],
    libdivide::divider<uint64_t> fast_prime,
    uint64_t p2 ,
    bool is_type_a,
    vector<uint64_t> &counts_backing_v,
    vector<uint64_t> &counts_backing_a,
    vector<uint64_t> &counts_backing_b,
    [[maybe_unused]] const long long int* cbv_negative_one,
    const long long int* cba_negative_one,
    const long long int* cbb_negative_one
)
{
    if (i_top < i_bot) return;

    uint64_t count = i_top - i_bot + 1;

    size_t j = 0;
    if (1) {
        // 99.8% in lower loops vs upper loop.
        if (count > 16) {
            __m256i v_ca = _mm256_set1_epi64x(c_a);
            __m256i v_cb = _mm256_set1_epi64x(c_b);

            for (; j+4 < count; j+=4) {
                counts[1] += 4;
                //size_t i = i_top - j;
                size_t i_low = i_top - j - 3;

                //const auto v = counts_backing_v[i];
                __m256i v_v = _mm256_loadu_si256((__m256i*)&counts_backing_v[i_low]);

                /* libdivide handles avx stuff */
                //uint64_t t = v / fast_prime;
                __m256i v_t = v_v / fast_prime;

                //size_t index = (t-1);
                __m256i v_index_plus_one = v_t;

#if INNER_ASSERT
                //assert(v >= p2);
                __m256i v_ge_p2 = _mm256_cmpgt_epi64(v_v, v_p2);
                assert(_mm256_movemask_epi8(v_ge_p2) == 0xFFFFFFFF);

                //assert(counts_backing_v[index] == t);
                __m256i v_t_check = _mm256_i64gather_epi64(cbv_negative_one, v_index_plus_one, 8);
                __m256i v_i_eq_t = _mm256_cmpeq_epi64(v_t_check, v_t);
                assert(_mm256_movemask_epi8(v_i_eq_t) == 0xFFFFFFFF);
#endif

                //uint64_t d_1 = counts_backing_a[index] - c_a;
                //uint64_t d_2 = counts_backing_b[index] - c_b;
                __m256i v_a_index = _mm256_i64gather_epi64(cba_negative_one, v_index_plus_one, 8);
                __m256i v_b_index = _mm256_i64gather_epi64(cbb_negative_one, v_index_plus_one, 8);
                __m256i v_d1 = _mm256_sub_epi64(v_a_index, v_ca);
                __m256i v_d2 = _mm256_sub_epi64(v_b_index, v_cb);

                //counts_backing_a[i] -= is_type_a ? d_1 : d_2;
                //counts_backing_b[i] -= is_type_a ? d_2 : d_1;
                __m256i v_d_a = is_type_a ? v_d1 : v_d2;
                __m256i v_d_b = is_type_a ? v_d2 : v_d1;
                __m256i v_a_i = _mm256_loadu_si256((__m256i*)&counts_backing_a[i_low]);
                __m256i v_b_i = _mm256_loadu_si256((__m256i*)&counts_backing_b[i_low]);
                __m256i v_res_1 = _mm256_sub_epi64(v_a_i, v_d_a);
                __m256i v_res_2 = _mm256_sub_epi64(v_b_i, v_d_b);

                _mm256_storeu_si256((__m256i*)&counts_backing_a[i_low], v_res_1);
                _mm256_storeu_si256((__m256i*)&counts_backing_b[i_low], v_res_2);
            }
        }
    }

    uint32_t i = i_top - j;
    for (; i >= i_bot; i--) {
        counts[2] += 1;
        const auto v = counts_backing_v[i];
        assert(v >= p2);

        uint64_t t = v / fast_prime;
        size_t index = (t-1);
        assert(counts_backing_v[index] == t);

        uint64_t d_1 = counts_backing_a[index] - c_a;
        uint64_t d_2 = counts_backing_b[index] - c_b;

        counts_backing_a[i] -= is_type_a ? d_1 : d_2;
        counts_backing_b[i] -= is_type_a ? d_2 : d_1;
    }
}

vector<uint64_t>
__get_special_prime_counts_vectorized_bulk(
        const uint64_t n, const uint32_t r,
        uint32_t start_prime,
        std::function< uint64_t(uint64_t)> init_count_a,
        std::function< uint64_t(uint64_t)> init_count_b,
        std::function< bool(uint64_t)> is_group_a
) {
    // Pair of how many numbers <= i of {form_a, form_b}
    // for i = 1, 2, ..., n/r, n/(r-1), n/(r-2), ... n/3 n/2 n/1

    /**
     * Could splitting counts_backing into [1, n/r-1] uint32_t and [n/r, n] uint64_t
     * Would save 25% memory
     * Indexing would be slightly easier and avx might be 2x faster (8 elements vs 4)
     */
    const auto length = ((n/r)-1) + r;
    vector<uint64_t> counts_backing_v;
    vector<uint64_t> counts_backing_a;
    vector<uint64_t> counts_backing_b;
    setup_vectors(n, r, init_count_a, init_count_b, counts_backing_v, counts_backing_a, counts_backing_b);
    assert(counts_backing_v.size() == length);

    // Do calculation | 98% of the work is here
    {
        primesieve::iterator it(/* start= */ start_prime);
        uint64_t prime = it.next_prime();
        assert(prime == start_prime);

        __m256i v_one = _mm256_set1_epi64x(1);

        /* These are shifted back 1 from the start of counts_backing_a.data()
         * so that we don't have to subtract 1 for the index */
        const long long int* cbv_negative_one = reinterpret_cast<const long long int*>(counts_backing_v.data()) - 1;
        const long long int* cba_negative_one = reinterpret_cast<const long long int*>(counts_backing_a.data()) - 1;
        const long long int* cbb_negative_one = reinterpret_cast<const long long int*>(counts_backing_b.data()) - 1;

        uint64_t counts[6] = {};
        const uint64_t K = 12;

        for (; prime <= r; prime = it.next_prime()) {
            libdivide::divider<uint64_t> fast_prime(prime);
            uint64_t p2 = prime * prime;
            __m256i v_p2 = _mm256_set1_epi64x(p2);

            const auto outer_v = counts_backing_v[prime-2];
            auto c_a = counts_backing_a[prime-2];
            auto c_b = counts_backing_b[prime-2];
            assert(outer_v == prime-1);

            // Index of last term < p2
            uint64_t stop_i = ((p2-1) < r) ? (p2-2) : (length - ((n-1) / p2 + 1));
            assert(counts_backing_v[stop_i] < p2);
            assert(counts_backing_v[stop_i+1] >= p2);

            /**
             * Updating counts of v in [p^2, n] V[k] -= (V[k/p] - V[p])
             */
            bool is_type_a = is_group_a(prime);

            /* i = (stop_i, length)
             * Update in reverse order
             *
             * First loop "Upper" loop handles when V[i] / prime > r
             * Second loop "Middle" loop uses avx
             *      Uses gather to find V[i] / prime
             * Third loop "Sequential" loop handles V[i] / prime = V[i+1] / prime = V[i+K] / prime
             *      V[i] / prime decreases sequentially hitting all numbers
             * Fourth loop "Bottom" handels i < r
             *      Exactly (V[i] / prime) is constant for p, update this many values with same value in a lopp
             */

            // N/j/prime >= r   <=>  N/j/prime >= N^(1/2)  <=>  j >= r/prime
            uint64_t i_threshold_1, i_threshold_2, i_threshold_3;
            {
                i_threshold_1 = length - n/(r*prime) - 1;
                assert( counts_backing_v[i_threshold_1] / prime < r );
                assert( counts_backing_v[i_threshold_1 + 1] / prime >= r );
                assert( length > i_threshold_1 );
                i_threshold_1 = std::max(i_threshold_1, stop_i);
                assert( i_threshold_1 >= stop_i );

                /* V[i] / prime - V[i+K] / prime < 1
                * n/j / prime - n/(j+K) / prime <= 1
                * n*(j+K) - n*(j) <= prime
                * n*K <= prime * j * (j+K)
                * j ~= sqrt(n * K / prime)
                */
                // Can disable by setting to r - 1;
                int64_t j = sqrt((double) n * K / prime);
                i_threshold_2 = (j < r) ? length - j : r;
                i_threshold_2 = std::max<uint64_t>(i_threshold_2, stop_i);

                i_threshold_3 = std::max<uint64_t>(stop_i, r - 1);

                //printf("p: %lu -> %lu, 1: %lu | 2: %lu | 3: %u | %lu\n", prime, length-1, i_threshold_1, i_threshold_2, r-1, stop_i);

                assert( i_threshold_1 >= i_threshold_2 );
                assert( i_threshold_2 >= stop_i );
                assert( i_threshold_2 >= r);
                assert( i_threshold_2 >= i_threshold_3);
            }


            // Upper loop handles i in (i_threshold_1, length) where V/prime > r
            for (size_t i = length - 1; i > i_threshold_1; i--) {
                counts[0] += 1;
                size_t index = length - (length - i) * prime;
                assert( counts_backing_v[i] / fast_prime == counts_backing_v[index] );

                uint64_t d_1 = counts_backing_a[index] - c_a;
                uint64_t d_2 = counts_backing_b[index] - c_b;

                counts_backing_a[i] -= is_type_a ? d_1 : d_2;
                counts_backing_b[i] -= is_type_a ? d_2 : d_1;
            }

            if (i_threshold_1 == stop_i) continue;


            // Middle loop handles i in (i_threshold_2, i_threshold_1], where i >= r
            middle_loop_avx(
                i_threshold_2 + 1, i_threshold_1,
                c_a, c_b,
                counts,
                fast_prime, p2, is_type_a,
                counts_backing_v, counts_backing_a, counts_backing_b,
                cbv_negative_one, cba_negative_one,  cbb_negative_one
            );

            // Sequential loop handles i in (i_threshold_3, i_threshold_2]
            {
                // dv_index = V[i] / prime - 1
                uint64_t dv = counts_backing_v[i_threshold_2] / fast_prime;
                // max value of V[i] dv is valid for
                uint64_t max_v = dv * prime;
                uint64_t index = dv - 1;

                uint64_t d_1 = counts_backing_a[index] - c_a;
                uint64_t d_2 = counts_backing_b[index] - c_b;
                uint64_t d_a = is_type_a ? d_1 : d_2;
                uint64_t d_b = is_type_a ? d_2 : d_1;

                for (size_t i = i_threshold_2; i > i_threshold_3; i -= 1) {
                    // Rarely happens, only once every K
                    if (counts_backing_v[i] < max_v) {
                        index -= 1;
                        max_v -= prime;
                        d_1 = counts_backing_a[index] - c_a;
                        d_2 = counts_backing_b[index] - c_b;
                        d_a = is_type_a ? d_1 : d_2;
                        d_b = is_type_a ? d_2 : d_1;
                    }

#if INNER_ASSERT
                    assert( counts_backing_v[i] / fast_prime == index + 1 );
#endif

                    counts_backing_a[i] -= d_a;
                    counts_backing_b[i] -= d_b;
                    counts[3] += 1;
                }
            }

            if (i_threshold_3 == stop_i) continue;

            // Bottom loop i < r
            if (1) {
                assert(i_threshold_3 == r-1);

                size_t i = r-1;

                // total number to process in this loop, stop when count == 0
                int64_t count;
                // v / prime, processing multiple at a time.
                size_t div;
                // Index of v[i] / prime
                size_t index;

                {
                    uint64_t v = r;
                    assert(counts_backing_v[i] == v);

                    count = i - stop_i;
                    div = v / fast_prime;
                    index = div - 1;
                }

                // Setup loop
                {
                    // v = r
                    int64_t mod = (int64_t) r - div * prime;
                    int64_t initial = count < mod ? count : mod;
                    //printf("p: %lu -> %ld @ %ld | start: %u, total: %ld\n", prime, initial, div, r, count);

                    assert(counts_backing_v[index] == div);
                    uint64_t d_1 = counts_backing_a[index] - c_a;
                    uint64_t d_2 = counts_backing_b[index] - c_b;

                    counts[4] += initial + 1;
                    for (; initial >= 0; i--, initial--) {
                        counts_backing_a[i] -= is_type_a ? d_1 : d_2;
                        counts_backing_b[i] -= is_type_a ? d_2 : d_1;
                    }
                    count -= initial;
                }

                // This loop is 28.1% of counts, <10% of time.
                for (; count >= prime; count -= prime) {
                    index--;

                    assert(counts_backing_v[i] / prime == counts_backing_v[index]);
                    uint64_t d_1 = counts_backing_a[index] - c_a;
                    uint64_t d_2 = counts_backing_b[index] - c_b;
                    uint64_t d_a = is_type_a ? d_1 : d_2;
                    uint64_t d_b = is_type_a ? d_2 : d_1;

                    for (size_t j = 0; j < prime; j++, i--) {
                        counts_backing_a[i] -= d_a;
                        counts_backing_b[i] -= d_b;
                    }
                    counts[5] += prime;
                }
                // Final loop with less than prime
                if (count > 0) {
                    index--;

                    assert(counts_backing_v[i] / prime == counts_backing_v[index]);
                    uint64_t d_1 = counts_backing_a[index] - c_a;
                    uint64_t d_2 = counts_backing_b[index] - c_b;

                    for (; count > 0; i--, count--) {
                        counts_backing_a[i] -= is_type_a ? d_1 : d_2;
                        counts_backing_b[i] -= is_type_a ? d_2 : d_1;
                    }
                    counts[4] += count;
                }
            } else {
                for (size_t i = i_threshold_3; i > stop_i; i--) {
                    counts[4] += 1;
                    const auto v = i - 1;
                    assert(v >= p2);
                    assert(v == counts_backing_v[i]);

                    uint64_t t = v / fast_prime;
                    size_t index = (t-1);
                    assert(counts_backing_v[index] == t);

                    uint64_t d_1 = counts_backing_a[index] - c_a;
                    uint64_t d_2 = counts_backing_b[index] - c_b;

                    counts_backing_a[i] -= is_type_a ? d_1 : d_2;
                    counts_backing_b[i] -= is_type_a ? d_2 : d_1;
                }
            }
        }
        fprintf(stderr, "\tLoop counts %lu, %lu, %lu, (K=%lu) %lu, %lu, %lu\n",
                counts[0], counts[1], counts[2], K, counts[3], counts[4], counts[5]);
    }

    return counts_backing_b;
}


vector<uint64_t>
__get_special_prime_counts_vectorized(
        const uint64_t n, const uint32_t r,
        uint32_t start_prime,
        std::function< uint64_t(uint64_t)> init_count_a,
        std::function< uint64_t(uint64_t)> init_count_b,
        std::function< bool(uint64_t)> is_group_a
) {
    // Pair of how many numbers <= i of {form_a, form_b}
    // for i = 1, 2, ..., n/r, n/(r-1), n/(r-2), ... n/3 n/2 n/1

    /**
     * Could splitting counts_backing into [1, n/r-1] uint32_t and [n/r, n] uint64_t
     * Would save 25% memory
     * Indexing would be slightly easier and avx might be 2x faster (8 elements vs 4)
     */
    const auto length = ((n/r)-1) + r;
    vector<uint64_t> counts_backing_v;
    vector<uint64_t> counts_backing_a;
    vector<uint64_t> counts_backing_b;
    setup_vectors(n, r, init_count_a, init_count_b, counts_backing_v, counts_backing_a, counts_backing_b);
    assert(counts_backing_v.size() == length);

    uint64_t counts[6] = {};

    // Do calculation | 98% of the work is here
    {
        primesieve::iterator it(/* start= */ start_prime);
        uint64_t prime = it.next_prime();
        assert(prime == start_prime);

        __m256i v_one = _mm256_set1_epi64x(1);

        /* These are shifted back 1 from the start of counts_backing_a.data()
         * so that we don't have to subtract 1 for the index */
        const long long int* cbv_negative_one = reinterpret_cast<const long long int*>(counts_backing_v.data()) - 1;
        const long long int* cba_negative_one = reinterpret_cast<const long long int*>(counts_backing_a.data()) - 1;
        const long long int* cbb_negative_one = reinterpret_cast<const long long int*>(counts_backing_b.data()) - 1;

        for (; prime <= r; prime = it.next_prime()) {
            libdivide::divider<uint64_t> fast_prime(prime);
            uint64_t p2 = prime * prime;
            __m256i v_p2 = _mm256_set1_epi64x(p2);

            const auto outer_v = counts_backing_v[prime-2];
            auto c_a = counts_backing_a[prime-2];
            auto c_b = counts_backing_b[prime-2];
            assert(outer_v == prime-1);

            // Index of last term < p2
            uint64_t stop_i = ((p2-1) < r) ? (p2-2) : (length - ((n-1) / p2 + 1));
            assert(counts_backing_v[stop_i] < p2);
            assert(counts_backing_v[stop_i+1] >= p2);

            /**
             * Updating counts of v in [p^2, n] V[k] -= (V[k/p] - V[p])
             */
            bool is_type_a = is_group_a(prime);

            /* N/j/prime >= r   <=>  N/j/prime >= N^(1/2)  <=>  j >= r/prime **/
            uint64_t i_break = length - n/(r*prime);
            assert( counts_backing_v[i_break-1] / prime < r );
            assert( counts_backing_v[i_break] / prime >= r );
            assert( i_break >= stop_i );

            if (0) {
                // Older, Slower
                for (size_t i = counts_backing_v.size() - 1; i > stop_i; i--) {
                    const auto v = counts_backing_v[i];
                    assert(v >= p2);

                    uint64_t t = v / fast_prime;
                    // size_t index = (t < r) ? (t-1) : (length - (n / t));
                    size_t index = (i < i_break) ? (t - 1) : (length - (n / t));
                    assert( (i < i_break) == (t < r) );
                    assert(counts_backing_v[index] == t);

                    if (0) {
                        // Older, Slower
                        if (is_type_a) {
                            counts_backing_a[i] -= counts_backing_a[index] - c_a;
                            counts_backing_b[i] -= counts_backing_b[index] - c_b;
                        } else {
                            counts_backing_a[i] -= counts_backing_b[index] - c_b;
                            counts_backing_b[i] -= counts_backing_a[index] - c_a;
                        }
                    } else {
                        uint64_t d_1 = counts_backing_a[index] - c_a;
                        uint64_t d_2 = counts_backing_b[index] - c_b;

                        counts_backing_a[i] -= is_type_a ? d_1 : d_2;
                        counts_backing_b[i] -= is_type_a ? d_2 : d_1;
                    }
                }
            } else {
                // Upper loop handes i >= i_break, where V/prime > r
                for (size_t i = counts_backing_v.size() - 1; i >= i_break; i--) {
                    size_t index = length - (length - i) * prime;
                    assert( counts_backing_v[i] / prime == counts_backing_v[index] );

                    uint64_t d_1 = counts_backing_a[index] - c_a;
                    uint64_t d_2 = counts_backing_b[index] - c_b;

                    counts_backing_a[i] -= is_type_a ? d_1 : d_2;
                    counts_backing_b[i] -= is_type_a ? d_2 : d_1;
                }
                uint64_t count = i_break <= stop_i ? 0 : i_break-1 - stop_i;
                middle_loop_avx(
                    stop_i + 1, i_break - 1,
                    c_a, c_b,
                    counts,
                    fast_prime, p2, is_type_a,
                    counts_backing_v, counts_backing_a, counts_backing_b,
                    cbv_negative_one, cba_negative_one,  cbb_negative_one
                );
            }
        }
    }

    return counts_backing_b;
}

vector<pair<uint64_t, pair<uint64_t, uint64_t>>>
__get_special_prime_counts(
        uint64_t n, uint32_t r,
        uint32_t start_prime,
        std::function< uint64_t(uint64_t)> init_count_a,
        std::function< uint64_t(uint64_t)> init_count_b,
        std::function< bool(uint64_t)> is_group_a
) {
    // Pair of how many numbers <= i of {form_a, form_b}
    // for i = 1, 2, ..., n/r, n/(r-1), n/(r-2), ... n/3 n/2 n/1

    vector<pair<uint64_t, pair<uint64_t, uint64_t>>> counts_backing;
    {
        size_t size = r + n/r - 1;
        counts_backing.reserve(size);
        // 1, 2, ... n / r - 1
        for(uint32_t v = 1; v < (n / r); v++) {
            uint64_t c_a = init_count_a(v);
            uint64_t c_b = init_count_b(v);
            assert(c_a + c_b <= v);
            counts_backing.push_back({v, {c_a, c_b}});
        }

        // n/r, n/(r-1), n/(r-2), ... n/3 n/2 n/1
        for(uint64_t i = r; i >= 1; i--) {
            uint64_t v = n / i;
            uint64_t c_a = init_count_a(v);
            uint64_t c_b = init_count_b(v);
            assert(c_a + c_b <= v);
            counts_backing.push_back({v, {c_a, c_b}});
        }
    }

    const auto length = ((n/r)-1) + r;
    assert(counts_backing.size() == length);

    {
        primesieve::iterator it(/* start= */ start_prime);
        uint64_t prime = it.next_prime();
        assert(prime == start_prime);
        for (; prime <= r; prime = it.next_prime()) {
            uint64_t p2 = prime * prime;

            auto& [v, c__] = counts_backing[prime-2];
            assert(v == prime-1);
            auto [c_a, c_b] = c__;

            // Index of last term < p2
            uint64_t stop_i = ((p2-1) < r) ? (p2-2) : (length - ((n-1) / p2 + 1));
            assert(counts_backing[stop_i].first < p2);
            assert(counts_backing[stop_i+1].first >= p2);

            bool is_type_a = is_group_a(prime);
            for (size_t i = counts_backing.size() - 1; i > stop_i; i--) {
                auto& [v, u] = counts_backing[i];
                assert(v >= p2);

                uint64_t t = v / prime;
                size_t index = (t < r) ? (t-1) : (length - (n / t));
                assert(counts_backing[index].first == t);

                const auto temp = counts_backing[index].second;
                if (0) {
                    if (is_type_a) {
                        u.first  -= temp.first  - c_a;
                        u.second -= temp.second - c_b;
                    } else {
                        u.first  -= temp.second  - c_b;
                        u.second -= temp.first   - c_a;
                    }
                } else {
                    uint64_t d_1 = temp.first - c_a;
                    uint64_t d_2 = temp.second - c_b;

                    u.first -= is_type_a ? d_1 : d_2;
                    u.second -= is_type_a ? d_2 : d_1;
                }
            }
        }
    }
    return counts_backing;
}


/**
 * Get number of primes <= i for important values of i.
 * Returns result in vector
 *
 * See get_special_prime_counts_map
 */
vector<uint64_t>
get_special_prime_counts_vector(
        uint64_t n, uint32_t r,
        uint32_t start_prime,
        std::function< uint64_t(uint64_t)> init_count_a,
        std::function< uint64_t(uint64_t)> init_count_b,
        std::function< bool(uint64_t)> is_group_a
) {
    auto start = std::chrono::high_resolution_clock::now();

    vector<uint64_t> count_primes;
    if (1) {
        if (1) {
            count_primes = __get_special_prime_counts_vectorized_bulk(
                n, r, start_prime,
                init_count_a, init_count_b, is_group_a);
        } else {
            count_primes = __get_special_prime_counts_vectorized(
                n, r, start_prime,
                init_count_a, init_count_b, is_group_a);
        }
    } else {
        // Older Slower
        auto counts_backing = __get_special_prime_counts(
            n, r, start_prime,
            init_count_a, init_count_b, is_group_a);

         /**
         * Reduces 3x uint64 to 1x uint64
         * BUT temporarily increases memory pressure
         *
         * Possibly I can read/erase from back, then reserve
         */
        count_primes.resize(counts_backing.size());
        for (size_t i = 0; i < counts_backing.size(); i++) {
            count_primes[i] = counts_backing[i].second.second;
        }
    }

    {
        auto end = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double>(end - start).count();
        fprintf(stderr, "\tcount_special_primes(%lu) = %lu  (%.2f)\n",
                n, count_primes.back(), elapsed);
    }

    // counts should be increasing
    assert(is_sorted(count_primes.begin(), count_primes.end()));

    return count_primes;
}

uint64_t isqrt(uint64_t x) {
    /* Power of 4 greater or equal to x */
    uint64_t q = 1;
    while (q <= x)
        q <<= 2;

    /* Quadaratic residue */
    uint64_t r = 0;
    while (q > 1) {
        q >>= 2;
        int64_t t = x - r - q;
        r >>= 1;
        if (t >= 0) {
            x = t;
            r += q;
        }
    }
    return r;
}

/**
 * Count population of 2^n of quadratic form generated by [2][p][q]^2
 *
 * TODO: can I calculate this directly using a variation of
 * get_special_prime_counts. Possibly computing "number less than n without
 * a duplicate factor".
 */
uint64_t count_population_quadratic_form(
        size_t bits,
        uint32_t start_prime,
        uint32_t add_to_special_primes,
        std::function< uint64_t(uint64_t)> init_count_a,
        std::function< uint64_t(uint64_t)> init_count_b,
        std::function< bool(uint64_t)> is_group_a
) {
    assert (bits <= 62); // overflows probably exist
                         //
    uint64_t n = 1ul << bits;
    uint64_t r = isqrt(n);
    assert(r*r <= n);
    assert((r+1) * (r+1) > n);
    assert(r < std::numeric_limits<uint32_t>::max());

    auto start = std::chrono::high_resolution_clock::now();

    // 50-65% of time is building special prime counts.
    const auto count_special_primes = get_special_prime_counts_vector(
        n, r,
        start_prime,
        init_count_a,
        init_count_b,
        is_group_a);

    const auto special_counts = count_special_primes.size();
    assert(special_counts == (n/r-1) + r);

    // TODO does (v <= r) work when n does not evenly divide r
    // TODO try inlining
    std::function<uint32_t(uint64_t)> count_special_index
      = [special_counts, n, r](uint64_t v) {
      return (v <= r) ? v - 1 : (special_counts - n / v);
    };

    // Only interested in these primes to odd powers
    vector<uint32_t> special_primes;
    {
        if (add_to_special_primes) {
          special_primes.push_back(add_to_special_primes);
        }

        size_t past = 0;
        primesieve::iterator it(/* start= */ start_prime);
        uint64_t prime = it.next_prime();
        assert(prime == start_prime);
        for (; past < 2; prime = it.next_prime()) {
            if (!is_group_a(prime)) {
                special_primes.push_back(prime);
                past += prime > r;
            }
        }
        assert(special_primes[special_primes.size() - 2] > r);  // Need two past r
        fprintf(stderr, "\tPrimes(%lu) = %u %u ... %u, %u, %u\n",
            special_primes.size(),
            special_primes[0],
            special_primes[1],
            special_primes[special_primes.size() - 3],
            special_primes[special_primes.size() - 2],
            special_primes[special_primes.size() - 1]);
    }

    std::function<uint64_t(uint64_t, uint32_t)> count_in_ex_large
      = [&special_primes, &count_special_primes, &count_special_index](uint64_t n, uint32_t pi) {
        // Generally happens when primes[pi-1] <= n < primes[pi]
        if (n < special_primes[pi]) {
            return n;
        }

        uint64_t count = n;

        // Handle p <= n < p^2
        uint64_t start_p = special_primes[pi];
        assert(start_p * start_p > n);

        uint32_t first_m = n / start_p;
        assert(first_m < start_p);

        uint64_t last = n / (first_m + 1);
        uint64_t count_last = 0;

        // Determine start prime of first loop of sqrt code.
        if (last < start_p) {
          assert(last <= special_primes.back());
          count_last = pi;
        } else {
          count_last = count_special_primes[count_special_index(last)];
        }

        // 75%+ of count_in_ex time is in this loop.
        // only 25% have first_m > 8
        for (uint32_t m = first_m; m > 0; m--) {
            // Count of number of primes with n / p == m
            //   -> Primes in the interval (n / (m + 1), n / m]
            // uint64_t first = last;
            last  = n / m;
            uint64_t count_first = count_last;

            count_last = count_special_primes[count_special_index(last)];

            assert(count_last >= count_first);
            count -= m * (count_last - count_first);
        }
        return count;
    };

    std::function<uint64_t(uint64_t, uint32_t, uint8_t)> count_in_ex;
    count_in_ex = [&special_primes, &count_in_ex, &count_in_ex_large]
                      (uint64_t n, uint32_t pi, uint8_t root) {
        if (n < special_primes[pi])
            return n;

        uint64_t count = 0;
        {
            /* Verify n < p^4, avoid overflow with p near uint32_t*/
            uint64_t p = special_primes[pi];
            uint64_t p2 = p * p;
            uint64_t n_p2 = n / p2;
            assert( (n_p2 < p2) || (n_p2 == p2 && (n % p2 > 0)));
        }

        assert(root <= 3);
        uint64_t fourth_root = isqrt(isqrt(n));

        if (root >= 3) {
          // Handle p where p^3 <= n < p^4
          // p^1 has recursion twice
          // p^2 has trivial recursion
          // p^3 has no recursion
          // Handle p^2 < n < p^3, only need to handle p^1 not p^3
          for (; pi < special_primes.size(); pi++) {
              uint64_t p = special_primes[pi];
              uint64_t p2 = p * p;

              uint64_t tn = n / p;
              if (p2 > tn)
                  break;

              // p^3 < n < p^4  ->  p^2 <= tn < p^3

              // TODO could possible pass a "skip first loop flag" (or call a
              // method that handles large p)
              count -= count_in_ex(tn, pi+1, 2);

              tn /= p; // p <= tn < p^2
              assert(tn >= p);

              // For each of these primes run the final loop logic
              count += count_in_ex_large(tn, pi+1);

              tn /= p; // 1 <= tn < p
              assert(tn < p);
              count -= tn;
          }
        }

        if (root >= 2) {
          // Handle p^2 <= n < p^3, only need to handle p^1 not p^3
          for (; pi < special_primes.size(); pi++) {
              uint64_t p = special_primes[pi];

              uint64_t tn = n / p;
              if (p > tn)
                  break;

              assert(tn >= p);
              count -= count_in_ex_large(tn, pi+1);
              // Have to add back all the counts of tn*r

              tn /= p;
              assert(tn < p);
              count += tn;  // count_in_exp(tn, pi+1);
          }
        }

        // Handles adding n to count.
        count += count_in_ex_large(n, pi);

        return count;
    };

    vector<pair<uint64_t,int32_t>> queue_n_pp_pi;
    /**
     * Handle when n >= p^4, e.g. root >= 4
     * Push to a vector which allows for parallelism in evaluating
     */
    std::function<uint64_t(int32_t sign, uint64_t, uint32_t)> count_in_ex_small;
    count_in_ex_small = [&special_primes, &queue_n_pp_pi, &count_in_ex_small]
                      (int32_t sign, uint64_t n, uint32_t pi) {
        // return doesn't use sign, queue does.
        if (n < special_primes[pi])
            return n;

        uint64_t count = 0;

        uint64_t fourth_root = isqrt(isqrt(n));

        // Handle p where p^4 <= n
        for (; pi < special_primes.size(); pi++) {
            uint64_t p = special_primes[pi];
            if (p > fourth_root) {
                queue_n_pp_pi.push_back({n, (int32_t) sign * pi});
                break;
            }

            uint64_t tn = n / p;

            // This loop has been optimized see A000047.py, for clearer code
            for (; ;) {
                if (tn < p) {
                    count -= tn;  // count_in_exp(tn, pi+1);
                    break;
                }

                count -= count_in_ex_small(-1 * sign, tn, pi+1);

                // Have to add back all the counts of tn * p
                tn /= p;
                if (tn < p) {
                    count += tn;  // count_in_exp(tn, pi+1);
                    break;
                }
                count += count_in_ex_small(sign, tn, pi+1);

                tn /= p;
            }
        }

        return count;
    };

    // in parallel handle all small primes
    uint64_t count = 0;
    if (1) {
        // Would be nice to break this into two sections with larger chunk in the later section.
        assert( special_primes.size() > 10 );
        auto start_1 = std::chrono::high_resolution_clock::now();

        // Break apart the small primes into individual powers, later handle larger groups at a time.
        count += count_in_ex_small(1, n, 0);
        // Sort largest n first so they get started first.
        sort(queue_n_pp_pi.rbegin(), queue_n_pp_pi.rend());
        fprintf(stderr, "\tbroke inclusion-exclusion into %lu cases\n", queue_n_pp_pi.size());

        auto start_2 = std::chrono::high_resolution_clock::now();

        #pragma omp parallel for reduction(+:count) schedule(dynamic, 2)
        for (const auto [n_pp, pi] : queue_n_pp_pi) {
            //printf("\tin_ex(%lu, %i(%u))\n", n_pp, pi, special_primes[abs(pi)]);
            if (pi >= 0)
                count += count_in_ex(n_pp, pi, 3);
            else
                count -= count_in_ex(n_pp, -pi, 3);
        }

        auto end = std::chrono::high_resolution_clock::now();

        double inex_split = std::chrono::duration<double>(start_2 - start_1).count();
        double inex_small = std::chrono::duration<double>(end - start_2).count();
        fprintf(stderr, "\tinclusion-exclusion took %.1f then %.1f\n", inex_split, inex_small);
    } else {
        count = count_in_ex(n, 0, 100);
    }

    {
        auto end = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double>(end - start).count();
        printf("| %2lu | %-18lu | %-18lu | %-8.2f |\n",
                bits, count, count_special_primes.back(), elapsed);
    }
    return count;
}
