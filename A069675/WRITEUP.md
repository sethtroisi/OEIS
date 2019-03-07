# A069675: "Primes with either no internal digits or all internal digits are 0.”

[TOC]

# Motivation

Recently good progress has been made on this sequence, In 2018 I discovered 100 (TODO EXACT COUNT) new terms.
These come from the combination of a O(log(sqrt(n))) = ~50x improvement to the sieving algorithm and increased access to compute.
I would appreciate some external eyes to validate my approach and suggestions on potential areas of improvement and investigation.

----

I use the shorthand (d,a,b) to represents the number a * 10^d + b, all candidates primes can be expressed in this form.

Note: I ignore all numbers less than ~1e11 (d <= 10) because it’s fast to test all 240 (d,a,b) tuples less than this limit and it’s nice for the candidate primes to be strictly larger than the primes I’m sieving with (otherwise terms in the small terms in the sequence may mark themselves off)

# 0. Helper Prework

Choose a MAX_D.
There are 24 * MAX_D candidates and they seem to be prime at the normal rate log(n = a * 10^d + b) ~= d * log(10) ~= d
I targeted 400 terms when I started and am now targeting 450 which requires MAX_D ~= 200,000

This number determines the number of expected primes and total cost of running. Both are estimated in [Real Numbers](#Real-Numbers)

# 1. Sieve

Let STATUS be a map of all candidates (d,a,b) to primes status
	in the code this is an array of size [MAX_D+1][10][10] initialized to -1 and then set to a prime factor if one is discovered.
	Two other special values are defined for “Prime” and “Composite with unknown factor”

Apply simple rules to mark candidates divisibility by 2,3,5, and 7

for (d,1,1)
  # Generalized Fermat Number
  test if d is a power of two if not mark candidate as divisible (with unknown number)
  # This could be slightly improved to mark off all numbers (as other people have tested all small values already)

This leaves 24 (a,b) pairs

    (1,1), (1,3), (1,7), (1,9), (2,3), (2,9), (3,1), (3,7), (4,3), (4,7), (4,9), (5,1),
    (5,3), (5,7), (6,1), (6,7), (7,1), (7,3), (7,9), (8,3), (8,5), (8,9), (9,1), (9,7)

Choose a large number, PRIME_CUTOFF, (more on choosing this later) and find all primes less using sieve of Eratosthenes. (in practice a segmented sieve is required)

    for each prime between 11 and PRIME_CUTOFF
        for all a,b pairs
            find a_inverse such that  a * a_inverse mod p = 1
            let test = a_inverse * -b % p
            observe that
                a * test + b mod p
                    = a * (a_inverse * - b) + b
                    = a * a_inverse * (-b)  + b
                    = 1 * -b + b
                    =  -b + b
                    = 0
            so a * test + b will be divisible by p

        let lookup be a map from all 24 values of test to their respective (a,b) pair
        the motivation for this is that now we can look at 10 ^ d mod p and if it matches any value in lookup then we know that it divides (d,a,b)

        for d > 10, d < MAX_D
            let t = 10 ^ d mod p
            if t in lookup:
                a,b from lookup[t]
                mark (d,a,b) divisible by p

This can be optimized further by noting that we do MAX_D iterations of 2nd loop which is >> 24

    choose value d_step, let it be close to sqrt(DMAX)
    find inverse_ten such that inverse_ten * 10 mod p = 1
    add an inner loop to calculate
        test * inverse_ten ^ d_i = (d_i, a, b)
    for 0 <= d_i <= d_step

modify the 2nd loop to increment by d_step and instead mark (d + d_i, a, b) divisible by p

stepup now takes O(24 * d_step) and the 2nd loop takes O(MAX_D / d_step)
this is minimized by setting d_step = sqrt(MAX_D / 24) and using a hash_table for lookup

	O(24 * d_increment * O(hash_inserts) + MAX_D / d_increment * O(hash_lookup))
	# hash_tables have theoretical O(1) insert and lookup
	# In practice on small values it more non constant but strongly sub-linear
	# Best real performance is found by setting d_increment a fraction (0.4-0.6) of it’s theoretical value

	minimize 24 * d_increment + MAX_D / d_increment
	# differentiate with respect to d_increment and set equal to 0
	0 = 24 - MAX_D / d_increment^2
	d_increment = sqrt(MAX_D / 24)

Theoretically this has reduced our work by a factor of sqrt(200,000/24) = 91
in practice this amounted to a 5-15x speedup

If MAX_PRIMES is chosen sufficiently large a great deal of the candidates primes are filtered.

    On a 32 thread Xeon machine with MAX_D = 200,000 the 4,800,000 candidates can be reduced to
	630,813 with MAX_PRIME=1e6 in 5 seconds
	551,595 with MAX_PRIME=1e7 in 40 seconds
    474,425 with MAX_PRIME=2e9 in 646 seconds


# 2. Testing

Given a list of candidate primes (d,a,b) perform Miller-Rabin primality test on each number
This is coded in C++ and leverages the high performance GMP library https://gmplib.org/manual/Prime-Testing-Algorithm.html
	Miller Rabin has no false negatives (marking a prime as composite) so any number marked composite is truly composite
	The chance of falsely accepting a number as prime is managed by setting REPS = 25
		The chance of adding a false prime is thus  < 1/4 ^ 25 * tests = 1e-15 * 1e6 = 1e-9
	Additionally the small list of ~300-500 primes can be verified with a much high REPS later.

Other than effecting managing partial status and resuming this code is trivially parallelizable and vastly more simple than sieving.
For MAX_D = 100,000, MAX_PRIME = 2e9 it took roughly 1000 cpu-core-days to test the 200,000 candidates
	This is roughly 7 minutes per number.

GMP makes use of FFT based multiplication so this code should have theoretic O(log^2(n)) = O(log^2(10^d)) = O(digits ^ 2)
In practice O(digits^2.3) better matches performance.

Hence the doubling of MAX_D from 100,000 to 200,000 will take roughly 10 times longer
	O(integral(d^2.2, d=0 to MAX_D)) = O(MAX_D ^ 3.3)
	This will likely yield 30 new terms and bring the total know primes to more than 450.
	Getting to 500 terms will require 100 the time and is sadly out of reach at this time.

# 3. Real Numbers

MAX_D = 200,000
Anecdotally 200 + log(D) primes => 450
![LogLog Plot](https://oeis.org/A069675/a069675_1.png)
A good prediction of count is 87 * log(log(10 * MAX_D)) =
So we expect 460 (MAX_D ~= 600,000 is needed for 500 terms)


PRIME_CUTOFF = 1e12
for my 32 thread Xeon sieve to this point this is expected to take ~4days and filter the 4,800,000 candidates down to ~350,000 candidates.
(TODO BETTER NUMBERS) Testing these number is expected to take 1XX days


# 4. Areas to investigate

* C++ performance of sieve has been extensively profiled.
  * Several different hash_sets, and hash_table approaches have been tested
    * The impact of these is large on the sieve step (2x performance) but ultimately only a 2-10% reduction in overall run time (as testing so vastly dominates)
  * Over 50%+ is spent in hash lookup
  * Further optimizations may exist but would need to be algorithmic (similar to recent d_step improvement)

* Trying to understand more general rules about divisibility of (a,b) pairs
  * Filtering (1,1) to d = power of two was a 5% improvement
  * Primes appear in all pairs but with potentially different frequency
    * There are 30 primes with a=7, b=9, 28 with a=6, b=1 but only 4 with a=4, b=7
    * This may be do to the number of candidates filtered by the sieve but a rule (similar to GFN) for any given (a,b) pair can be up a substantial improvement

* Out of curiosity I tried to move the sieve to a GPU with CUDA but the large amount of modules meant it was MUCH slower than the CPU version.
  * This seems in contrast to other results (e.g. https://github.com/curtisseizert/CUDASieve) and might be worth a 2nd look


# Appendix

Counts of Prime (a,b) pairs

(7, 9): 30, (9, 7): 28, (6, 1): 27, (8, 9): 25, (1, 3): 22, (1, 7): 22, (5, 3): 21, (7, 1): 21,
(6, 7): 21, (1, 9): 20, (9, 1): 19, (2, 3): 19, (5, 9): 18, (7, 3): 18, (3, 7): 17, (4, 3): 17,
(4, 1): 16, (3, 1): 16, (4, 9): 13, (2, 9): 8, (8, 3): 6, (4, 7): 4, (1, 1): 1

