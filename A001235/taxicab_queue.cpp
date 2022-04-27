#include<cassert>
#include<chrono>
#include<cstdio>
#include<queue>

struct item {
    long a3b3;
    long b;
};


class CompareItem {
    public:
    bool operator()(item &i1, item &i2)
    {
        return i1.a3b3 > i2.a3b3;
    }
};

int main(void)
{
    auto start = std::chrono::steady_clock::now();

    int  count = 0;
    long iters = 0;

    /**
     * This program is majorily limitted by O(log(n)) inserts in items
     * Extra branches / logic doesn't slow us down much
     *
     * On a Ryzen 3900x -> 12 - 15 million iters/second (pop + push)
     * 50,000 th pair (1601017842872) ~3.5 seconds
     * 1,000,000 th pair (1376426637022528) ~400 seconds
     */
    std::priority_queue<item, std::vector<item>, CompareItem> items;

    item i;
    items.push({1, 0}); // (1, 0) -> 1

    long next_a = 2;
    long next_a3 = next_a * next_a * next_a;

    long lastlast = -1;
    long last = -1;
    while (++iters)
    {
        i = items.top();
        items.pop();

        // Need to insert (next_a, 0)
        if (i.b == 0)
        {
            items.push({next_a3, 0});
            next_a++;
            next_a3 = next_a * next_a * next_a;
        }

        if (i.a3b3 == last && last != lastlast)
        {
            count++;
            if ((count <= 10)
                    || (count <= 1500 && count % 100 == 0)
                    || (count <= 50000 && count % 2000 == 0)
                    || (count % 10000 == 0)) {
                auto end = std::chrono::steady_clock::now();
                double elapsed = std::chrono::duration<double>(end-start).count();
                printf("  %7dth %17lu  (b: %4lu size: %lu  time: %.1f  iters: %lu)\n",
                        count, i.a3b3, i.b, items.size(), elapsed, iters);
            }

            if (count == 1000000)
                break;
        }
        lastlast = last;
        last = i.a3b3;

        // (a, b) -> (a, b+1)
        // (b+1)^3 - b^3 = 1 + 3*b + 3*b^2 = 1 + 3*(b + b*b)
        long b2 = i.b * i.b;
        long b3 = i.b * b2;
        if ( (b3 << 1) < i.a3b3 )
        {
            long b3_delta = 1 + 3ul * (i.b + b2);
            i.b++;
            i.a3b3 += b3_delta;
            items.push(i);
        }

    }
    return 0;
}
