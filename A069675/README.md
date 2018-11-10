[A069675](https://oeis.org/A069675)
---

### Background

I initially discovered A069675 as part of some Project Euler and took up the task of increasing the known bounds. In 2016 time Charles R Greathouse IV and Robert Israel had calculated up to a(300). I've extended this to a(434) and hope to eventually find a(450) which is likely to be around 10<sup>150,000</sup>

The code runs in two steps

1. Sieve all tuples `(a, b, d) = a * 10^d + b` with primes up to X million (10,000M for d <= 100,000).
1. gmp is_prime on all remaining pairs in parallel with checkpoints.


The largest run I've completed is `d <= 100,000` with an initial sieve of `primes <= 10B`.

```
git clone https://github.com/sethtroisi/OEIS
cd OEIS/A069675
git submodule update --init --recursive

sudo apt-get update
sudo apt-get install libgmp-dev

mkdir -p build
cd build

cmake ..
cmake --build . -- -j4

time ./sieve
time ./tester
```

`python A069675_helper.py` provides estimates of a(n), a guess at the optimal sieve size, and estimated wall time.

TODO table of results

![LogLog plot](results/LogLogA069675.png)

