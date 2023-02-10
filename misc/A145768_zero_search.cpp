#include <iostream>
using namespace std;

int main() {
    __int128_t x = 0;
    __int128_t i2 = 0;
    for (uint64_t i = 1; i < 100'000'000'000ul; i += 4) {
        // odd, odd, even, even
        // only need to look at  two even values
        i2 += 2*i - 1;
        x ^= i2;

        i2 += 2*i + 1;
        x ^= i2;

        i2 += 2*i + 3;
        x ^= i2;
        if (x == 0)
            cout << i+2 << endl;

        i2 += 2*i + 5;
        x ^= i2;
        if (x == 0)
            cout << i+3 << endl;
    }
}
