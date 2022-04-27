#include <cstdint>
#include <chrono>
#include <cstdio>
#include <iostream>
#include <queue>
#include <vector>

using std::priority_queue;
using std::vector;

using std::cout;
using std::endl;

struct data
{
    uint64_t n_3x2p4y2;
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

/** Population of 3 x^2 + 4 y^2 */
int main(void)
{
    auto start = std::chrono::steady_clock::now();

    priority_queue<data, std::vector<data>, compGT> items;

    uint64_t last_n = 0xFFFFFFFF; // so (0,0) doesn't match
    size_t count = 0;
    size_t iters = 0;

    uint64_t next_x, next_3x2;
    data item;

    {
        item.y = 0;
        item.n_3x2p4y2 = 0; // (0, 0)
        items.push(item);
    }

    next_x = 1;
    next_3x2 = 3 * next_x * next_x;

    uint32_t bits = 4;
    uint64_t next_log = 1ull << bits;

    long long operations = 0;

    while (true) //(!items.empty())
    {
        item = items.top();
        if (item.n_3x2p4y2 > next_3x2)
        {
            // replace temp item with (next_x, 0); the other item stays on the queue
            item.y = 0;
            item.n_3x2p4y2 = next_3x2;

            next_x++;
            next_3x2 = 3ull * next_x * next_x;
        }
        else
            items.pop();

        if (item.n_3x2p4y2 > next_log)
        {
            auto end = std::chrono::steady_clock::now();
            double elapsed = std::chrono::duration<double>(end-start).count();
            // Subtract 1 for 0
            printf("| %2d | %-11lu | %-13lu | %-7.2f secs | size: %5lu, iters/s: %.1f million \n",
                    bits, count - 1, iters, elapsed, items.size(), iters / 1e6 / elapsed);
            //if (bits == 33)
            //    break;
            bits += 1;
            next_log = 1ull << bits;
        }

        if (item.n_3x2p4y2 == last_n)
        {
            count--;

            // Increment all items with same n
            while (items.top().n_3x2p4y2 == last_n)
            {
                iters++;
                data tempItem = items.top(); items.pop();
                tempItem.n_3x2p4y2 += ((tempItem.y << 1) + 1) << 2;
                tempItem.y++;
                items.push(tempItem);
            }
        }

        last_n = item.n_3x2p4y2;
        item.n_3x2p4y2 += ((item.y << 1) + 1) << 2;
        item.y++;
        count++;
        iters++;
        items.push(item);
    }
}
