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
cat t | awk '{$16=$18=""; t=$15; $15=$17; $17=t; print $0}' | sed 's/  \+/ /g'
```


## Improvements

1. Lehmer uses `Fib(K)` as a lower bound on the size of the continued fraction.
   Something like `exp(sum(log(k_i)))` would enable avoiding 75%+ of expansion.

2. Continued fractions are recomputed for each subsequent `p`
   * You can save the solutions to Pell Equations but that requires expanding all
     of the continued fractions and not just using the Fib(K) trick

