# [A217259](https://oeis.org/A217259)

Should maybe be moved under [A050507](https://oeis.org/A050507)

This short line laid down a challange that cost me two nights of sleep.

> No term found below 2 \* 10^9 to continue sequence 435, 8576, 8826, ... - Michel Marcus, Mar 19 2013

code has now been tested up to `2e13 = 20,000,000,000,000`

## Double Check

All composite terms are verified against the know sequences:

[A050507](https://oeis.org/A050507) for dist=2
[A054903](https://oeis.org/A054903) for dist=3
[A063680](https://oeis.org/A063680) for dist=7
[A059118](https://oeis.org/A059118) for dist=8
[A054902](https://oeis.org/A054902) for dist=12

`primesieve` provides a near instant verification of prime count and dist=2 count.

`verify` can verify all the prime distance counts.


### Misc

Because `Sigma(i) = i + 1` iff i is prime, `is_prime(i) = (Sigma(i) - i == 1)`

Can use `sigma[i]` to check primality.

### Odd distances

There are only a handful of terms with odd distance.

* DIST=7 ([A063680](https://oeis.org/A063680)) with terms related to `2*3^n` and `2*3^n-7 = 2 * prime`
* These distances have only a single known term (and it's small)
  * DIST=13: single term (4418, 4431),  `4418=2*47^2`
  * DIST=17: single term (16, 33), `16=4^2`
  * DIST=19: single term (6, 25), `25=5^2`
  * DIST=37: single term (44, 81), `81=9^2`

In [A015886](https://oeis.org/A015886) T. D. Noe notes

> The "other" values of n are the odd n such that n+2 is not prime. For these n, in order for sigma(k) or sigma(n+k) to be odd, either k or n+k must be a square or twice a square. Examples: for n=7, n+k=9^2; for n=13, k=2\*47^2 and for n=19, n+k=5^2. Using this idea, it is easy to show that if a(23) exists it is greater than 10^12. - T. D. Noe, Sep 24 2007

This insight comes from inspecting the product definition of Sigma.

`1 + p + p^2 + ... + p^e = (p^(e+1) + 1) / (p-1)`

This is odd for all power of two and even powers of all other primes.
An odd `sigma(k)` implies that `k = 2^{0,1} * a^2` "Squares and twice squares" [A028982](https://oeis.org/A028982).

From the sequence definitian `sigma(k + odd) == sigma(k) + odd`. So one of `sigma(k + odd)` or `sigma(k)` is odd.
So all we need to look at is terms near squares (and twice squares).

Have to look both above squares, `(16, 16+17)`, and below, `(6+19, 25)`.

Have search up to 10^16 with `pypy check_squares.py` which only takes 50minutes
