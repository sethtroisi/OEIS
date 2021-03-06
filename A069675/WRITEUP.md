# [A069675](https://oeis.org/A069675) "Primes where all internal digits are 0” (e.g. 10007)

## Motivation

Recently good progress has been made on this sequence.

In 2018 I computed 130 new terms by sieving the candidates and distribution the computation on a large Xeon machine.

In 2019 I improved the sieving algorithm by factor of log(n) approximately a 50x improvement

#### Note

I use the shorthand `(d,a,b)` to represents the number `a * 10^d + b`, all candidates primes can be expressed in this form. I try to refer to this as a candidate and occasionally as a `(d,a,b)` tuble.

I ignore all numbers less than `~1e20 (d <= 20)` because it’s fast to test all small `(d,a,b)` candidates and it’s nice for the candidate primes to be strictly larger than the primes used for sieving with.

## 0. Helper Prework

#### Choose a `MAX_D`.

There are `24 * MAX_D` candidates and they seem to be prime at the normal rate `log(n = a * 10^d + b) ~= d * log(10) ~= d`<br>
I targeted 400 terms when I started and am now targeting 450 which requires `MAX_D = 200,000`

This number determines the number of expected primes and total cost of running. Both are estimated in [Real Numbers](#3-real-numbers)

## 1. Sieve

Let STATUS be a map of all candidates `(d,a,b)` to primes status (prime, composite, known prime factor)<br>

Apply simple rules to mark candidates divisibility by 2,3,5, and 7

For candidates of the form `(d,1,1)` (Generalized Fermat Number) test if d is a power of two if not mark candidate as divisible (composite)
(This could be slightly improved to mark off all numbers (as other people have tested all small values already) but it would only save testing ~20 candidates)

This leaves 22 (a,b) pairs

    (1,3), (1,7), (1,9), (2,3), (2,9), (3,1), (3,7), (4,3), (4,7), (4,9), (5,3),
    (5,9), (6,1), (6,7), (7,1), (7,3), (7,9), (8,3), (8,5), (8,9), (9,1), (9,7)

Choose a large number, `PRIME_CUTOFF`, (more on how to choosing this later) and find all primes less using sieve of Eratosthenes. (in practice a segmented sieve is required)

    for each prime between 11 and PRIME_CUTOFF
        for all a,b pairs
            find a_inverse such that a * a_inverse mod p = 1
            let test = (a_inverse * -b) mod p
            observe that
                a * test + b mod p
                    = a * (a_inverse * -b) + b
                    = a * a_inverse * (-b)  + b
                    = 1 * -b + b
                    = 0
            so a * test + b will be divisible by p

        let lookup be a map from all 24 values of test to their respective (a,b) pair
        the motivation for this is that now we can look at 10 ^ d mod p and if it matches any value in lookup then we know that p divides (d,a,b)

        for 10 < d <= MAX_D
            let test = 10 ^ d mod p
            if t in lookup:
                a,b from lookup[test]
                mark (d,a,b) divisible by p

This can be optimized further by noting that we do `MAX_D` iterations of 2nd loop

 choose a value `d_step`, let it be close to `sqrt(MAX_D)`

    find inverse_ten such that inverse_ten * 10 mod p = 1
    for 0 <= d_i <= d_step
        lookup[(d_i, a, b)] = test * inverse_ten ^ d_i

in the original code modify the 2nd loop to increment by `d_step` and instead mark `(d + d_i, a, b)` divisible by p

stepup now takes `O(24 * d_step)` and the 2nd loop takes `O(MAX_D / d_step)`
this is minimized by setting `d_step = sqrt(MAX_D / 24)` and using a hash_table for lookup

    O(24 * d_step * O(hash_inserts) + MAX_D / d_step * O(hash_lookup))

hash tables have theoretical `O(1)` insert and lookup, in practice on small values it's more complicated.

minimize `24 * d_step + MAX_D / d_step` by differentiate with respect to `d_step` and setting equal to 0

    0 = 24 - MAX_D / d_step^2
    d_step = sqrt(MAX_D / 24)

Best real world performance is found by setting `d_step` to half of it’s theoretical value.

Theoretically this has reduced our work by a factor of `sqrt(200,000/24) = 91`<br>
in practice this amounted to a 10-50x speedup.

If `MAX_PRIMES` is chosen sufficiently large a great deal of the candidates primes are filtered. On a 32 thread Xeon machine with `MAX_D = 200,000` the 4,800,000 candidates can be reduced to

    630,813 with MAX_PRIME=1e6 in 5 seconds
    551,595 with MAX_PRIME=1e7 in 40 seconds
    474,425 with MAX_PRIME=2e9 in 646 seconds

## 2. Testing

Given a list of candidates `(d,a,b)` not filtered by the sieve perform Miller-Rabin primality test on each number. This is coded in C++ and leverages the highly performant [GMP library](https://gmplib.org/manual/Prime-Testing-Algorithm.html).

Miller Rabin has no false negatives (marking a prime as composite) so any number marked composite is truly composite. The chance of falsely accepting a number as prime is managed by setting `REPS = 25`, the chance of adding a false prime is thus  `< 1/4 ^ REPS * number_of_candidates = 1/4^25 * 1e6 = 1e-9`.

Other than effectively managing checkpoints and resuming this code is trivially parallelizable and vastly more simple than the sieving code.
For `MAX_D = 100,000, MAX_PRIME = 2e9` it took roughly 1000 cpu-core-days to test the ~200,000 candidates (roughly 7 minutes per number).

GMP makes use of FFT based multiplication so this code should have theoretic `O(log^2(n)) = O(log^2(10^d)) = O(digits ^ 2)` In practice `O(digits^2.3)` better matches performance.

Hence the doubling of `MAX_D` from 100,000 to 200,000 will take roughly 10 times longer

    O(integral(d^2.3, d=0 to MAX_D)) = O(MAX_D ^ 3.3)

This will likely yield 30 new terms and bring the total know primes to more than 450.
Getting to 500 terms will require 100 the time and is sadly out of reach at this time.

## 3. Real Numbers

A great prediction of count is `36 * HarmonicNumber(X)`<br>
This is equal to `21 + 36 * log(MAX_D)` for reasonable values of `MAX_D`

So we expect 460 from `MAX_D = 200,000` and 500 from `MAX_D = 600,000`
![LogLog Plot](https://oeis.org/A069675/a069675_1.png)

`PRIME_CUTOFF = 1e12` for my 32 thread Xeon his is expected to take 4 days and filter the 4,800,000 candidates down to ~370,000 candidates.

Testing these number is expected to take between 60 and 200 days (better estimate when testing starts)

----

One interesting note is that sieve removes ~50% more candidates (than a small brute force sieve)<br>

    integral(d^2.3, d=0 to MAX_D) = integral(1/2 * d^2.3, d=0 to X*MAX_D)
    MAX_D^3.3 = 1/2 * (X*MAX_D)^3.3
    MAX_D^3.3 = 1/2 * X^3.3 * MAX_D^3.3
    X^3.3 = 2
    X = 2 ^ -3.3 = 1.234

This means we can test a 23% higher `MAX_D`, which yields on average 5.8 more primes.

## 4. Areas to investigate

* C++ performance of sieve has been extensively profiled.
  * Several different hash sets and hash table approaches have been tested
     * The impact of these is large on the sieve step (2x performance) but only a 2-10% reduction in overall run time (as testing so vastly dominates)
  * Over 50%+ is spent in hash lookup
  * Further optimizations may exist but would need to be algorithmic (similar to recent `d_step` improvement)

* Trying to understand more general rules about divisibility of (a,b) pairs
  * Filtering (1,1) to d = power of two was a 5% improvement
  * Primes appear in all pairs but with potentially different frequency
    * There are 30 primes of the form `(d,7,9)`, 28 of the form `(d,6,1)` but only 4 of the form `(d,4,7)`
    * This may be do to the number of candidates filtered by the sieve but a rule (similar to GFN) for any given (a,b) pair can be up a substantial improvement

* Out of curiosity I tried to move the sieve to a GPU with CUDA but the large amount of modules meant it was MUCH slower than the CPU version.
  * This seems in contrast to other results (e.g. https://github.com/curtisseizert/CUDASieve) and might be worth a 2nd look


# Appendix

Counts of Prime (a,b) pairs

    (7, 9): 30, (9, 7): 28, (6, 1): 27, (8, 9): 25, (1, 3): 22, (1, 7): 22, (5, 3): 21, (7, 1): 21,
    (6, 7): 21, (1, 9): 20, (9, 1): 19, (2, 3): 19, (5, 9): 18, (7, 3): 18, (3, 7): 17, (4, 3): 17,
    (4, 1): 16, (3, 1): 16, (4, 9): 13, (2, 9): 8, (8, 3): 6, (4, 7): 4, (1, 1): 1

Counts of (a,b) pairs that weren't filtered

    (7,9) : 22553
    (9,7) : 22503
    (7,1) : 21550
    (1,9) : 21376
    (9,1) : 21321
    (1,7) : 21318
    (8,9) : 19882
    (5,9) : 19806
    (6,7) : 19487
    (5,3) : 18211
    (6,1) : 18204
    (2,3) : 17296
    (4,9) : 16041
    (3,1) : 15315
    (1,3) : 15288
    (3,7) : 14592
    (7,3) : 14335
    (4,3) : 13908
    (4,1) : 11817
    (8,3) : 10245
    (4,7) : 6758
    (2,9) : 6403
    (1,1) : 15

Relative rate of primes (primes / tested) ignoring (1,1)

    (6, 1): 27/18204 prime = 0.0015
    (4, 1): 16/11817 prime = 0.0014
    (1, 3): 22/15288 prime = 0.0014
    (8, 9): 25/19882 prime = 0.0013
    (7, 9): 30/22553 prime = 0.0013
    (7, 3): 18/14335 prime = 0.0013
    (9, 7): 28/22503 prime = 0.0012
    (5, 3): 21/18211 prime = 0.0012
    (4, 3): 17/13908 prime = 0.0012
    (3, 7): 17/14592 prime = 0.0012
    (2, 9): 8/6403   prime = 0.0012 (d = 5, 25, 455, 761, 9205, 13561, 15955, 26669)
    (6, 7): 21/19487 prime = 0.0011
    (2, 3): 19/17296 prime = 0.0011
    (7, 1): 21/21550 prime = 0.0010
    (3, 1): 16/15315 prime = 0.0010
    (1, 7): 22/21318 prime = 0.0010
    (9, 1): 19/21321 prime = 0.0009
    (5, 9): 18/19806 prime = 0.0009
    (1, 9): 20/21376 prime = 0.0009
    (4, 9): 13/16041 prime = 0.0008
    (8, 3): 6/10245  prime = 0.0006 (d = 31, 105, 113, 369, 1359, 6219)
    (4, 7): 4/6758   prime = 0.0006 (d = 3, 9, 39, 2323)

Sieving might be improved for rules in "Primes of the Form" (TODO link)?

    # http://people.math.umass.edu/~bates/Primes_of_the_form.pdf

    # The most famous result
    p = x^2 + y^2 <=> p congruent to 1 mod 4

    p = x^2 + 2y^2 <=> p congruent to 1,3 mod 8
    p = x^2 + 3y^2 <=> p congruent to 1,7 mod 12
    p = x^2 + 4y^2 <=> p congruent to 1 mod 4
    p = x^2 + 7y^2 <=> p congruent to 1,9,11,15,23,25 mod 4

    p = x^2 + 6y^2 <=> p  congruent to 1,7 mod 24
    p = x^2 + 10y^2 <=> p congruent to 1,9,11,19 mod 40
    p = x^2 + 30y^2 <=> p congruent to 1,31,49,79 mod 120

    # These can be mapped to rules about a candidate (d,a,b)
    # e.g. (2*c, 4, 1)
    #       = 4 * 10 ** (2**c) + 1
    #       = 4 * (10**c)^2 + 1^2
    #   if this is prime then it will be 1 mod 4
    #   if (2*c, 4, 1) % 4 was not congruent to 1 the canditate could be
    #   be filtered immediatly

    #   So this didn't help :(
