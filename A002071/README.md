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

1. PQa to solve Pell Equations is 2-3x slower.
   * The continued fraction (CF) approach gets to do most of it's math in `__uint128_t`
      * The `uint128_t` code is 10x faster than the `mpz_class` code!
      * At `P=73` with `log2(prod(primes(73))) = 95`, 90% of `D` are < 64 bits.
      * At `P=97` with `log2(prod(primes(73))) = 95`, 50% of `D` are < 64 bits.
      * At `P=131` >99.7% of `D` are < 127 bits and on the fast path!
      * Maybe at `P=151` we'd have a few percent but that's still a LOT of fast path.

1. Found should be written to disk. Maybe even working Pell solutions
   * This goes back to the idea of expanding not redoing work for larger P OR of doing all the passes in 1 pass.
