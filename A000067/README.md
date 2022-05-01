# [A00067](https://oeis.org/A00067)

This sequence is trivially similiar (in solution) to [A000047](../A000047)/[A000050](../A000050)

## Results

TBD

| n  | a(n)          | count p in (5,7) mod 8 <= 2^n | time (s) |

```
# Fast
$ g++ -O3 -std=c++17 A000067_signature.cpp -lprimesieve

# For writing out the b-file
$ cat README.md  | awk '/\| .. \| ./ {print $2, $4}'
```
