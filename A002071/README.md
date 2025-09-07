# [A002071](https://oeis.org/A002071)

Related
[A002072](https://oeis.org/A002072)
and
[Project Euler 581](https://projecteuler.net/problem=581)

[Blog Post with larger records](https://11011110.github.io/blog/2007/03/23/smooth-pairs.html)
And
[Lehmer's "On a Problem of Stormer's" Paper](https://scispace.com/pdf/on-a-problem-of-stormer-2mpsudvm07.pdf)

## Compile and run

```
g++ -O2 --std=c++23 -o A002071 A002071.cpp -l gmpxx -lgmp -fopenmp
time ./A002071 43
```

## Comparison

```
cat t | awk '{$16=$18=""; print $0}' | sed 's/  \+/ /g'
```


## Improvements

1. Lehmer uses `Fib(K)` as a lower bound on the size of the continued fraction.
   Something like `exp(sum(log(k_i)))` would enable avoiding 75%+ of expansion.

2. Continued fractions are recomputed for each subsequent `p`
   * You can save the solutions to Pell Equations but that requires expanding all
     of the continued fractions and not just using the `Fib(K)` trick
   * Given that we know what P we are running up to (and hence `Max(MAX_CF)`) we
     could store `(D, CF (list))` for the X million `D` that generate `CF` that could
     possibly be smooth.
     * Not clear but maybe `x_1` is smaller than `CF (list)`? also there are so few `x_1`
       that you can just check if they are max(P)-smooth and only store them if so. This gives
       both fast partial answers for each `Q` but also avoids recomputing for each `P`!
   * Basically always do the expansion for `max(P)` but then save results for later (or
       just mark them into `Allstats[P]`).

3. PQa to solve Pell Equations is 2-3x slower.
   * The continued fraction (CF) approach gets to do most of it's math in `__uint128_t`
      * The `uint128_t` code is 10x faster than the `mpz_class` code!
      * At `P=73` with `log2(prod(primes(73))) = 95`, 90% of `D` are < 64 bits.
      * At `P=97` with `log2(prod(primes(73))) = 95`, 50% of `D` are < 64 bits.
      * At `P=131` >99.7% of `D` are < 127 bits and on the fast path!
      * Maybe at `P=151` we'd have a few percent but that's still a LOT of fast path.

4. Found should be written to disk. Maybe even working Pell solutions
   * This goes back to the idea of expanding not redoing work for larger P OR of doing all the passes in 1 pass.
