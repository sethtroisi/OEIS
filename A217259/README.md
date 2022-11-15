# [A217259](https://oeis.org/A217259)

Should maybe be moved under [A050507](https://oeis.org/A050507)

This short line laid down a challange that cost me two nights of sleep.

> No term found below 2 \* 10^9 to continue sequence 435, 8576, 8826, ... - Michel Marcus, Mar 19 2013

All even distances (up to 10) have now been tested up to `25.6e12 = 25,600,000,000,000`

All odd distances (up to 19) have been tested up to `1e17 = 100,000,000,000,000,000`

## Double Check

All composite terms are verified against the know sequences:

[A050507](https://oeis.org/A050507) for dist=2
[A054903](https://oeis.org/A054903) for dist=3
[A063680](https://oeis.org/A063680) for dist=7
[A059118](https://oeis.org/A059118) for dist=8
[A054902](https://oeis.org/A054902) for dist=12

`primesieve` provides a near instant verification of prime count and dist=2 count.

`verify.py` can verify all the prime distance counts.


### Misc

Because `Sigma(i) = i + 1` iff i is prime, `is_prime(i) = (Sigma(i) - i == 1)`

Can use `sigma[i]` to check primality.

### Odd distances

There are only a handful of terms with odd distance.

* DIST=7 ([A063680](https://oeis.org/A063680))
  * terms are related to `2*3^n` and `2*3^n-7 = 2 * prime`
  * **17 known terms!**
* DIST=1097 with **two known terms**
  * (8836, 9933), `8836=94^2`
  * (227529, 228626), `227529=477^2`
* These are small distances (all have only a single known term, despite searches up to 10^17)
  * DIST=13: (4418, 4431),  `4418=2*47^2`
  * DIST=17: (16, 33), `16=4^2`
  * DIST=19: (6, 25), `25=5^2`
  * DIST=37: (44, 81), `81=9^2`
  * DIST=39: (10, 49), `49=7^2`
  * DIST=49: (3315, 3364), `3364=58^2`
  * DIST=77: (98, 175), `98=2*7^2`
  * DIST=93: (72, 165), `72=2*6^2`
* high water marks for square terms (excluding DIST=7)
  * DIST=13:	(4418, 4431), `4418=2*47^2`
  * DIST=49:	(3315, 3364), `3364=58^2`
  * DIST=161:	(1435204, 1435365), `1435204=1198^2`
  * DIST=831:	(9017178, 9018009), `9018009=3003^2`
  * DIST=8095:	(435590546, 435598641), `435598641=20871^2`
  * DIST=77849:	(3078696196, 3078774045), `3078696196=55486^2`
  * DIST=98525:	(3341649249, 3341747774), `3341649249=57807^2`

* All distances with known terms
  * 7, 13, 17, 19, 37, 39, 49, 77, 93
  * 137, 161, 217, 221, 237, 245, 327, 329, 413, 421, 427,
  * 529, 573, 649, 655, 739, 757, 831, 889, 903, 937, 957, 977

In [A015886](https://oeis.org/A015886) T. D. Noe notes

> The "other" values of n are the odd n such that n+2 is not prime. For these n, in order for sigma(k) or sigma(n+k) to be odd, either k or n+k must be a square or twice a square. Examples: for n=7, n+k=9^2; for n=13, k=2\*47^2 and for n=19, n+k=5^2. Using this idea, it is easy to show that if a(23) exists it is greater than 10^12. - T. D. Noe, Sep 24 2007

This insight comes from inspecting the product definition of Sigma.

`1 + p + p^2 + ... + p^e = (p^(e+1) + 1) / (p-1)`

This is odd for all power of two and even powers of all other primes.
An odd `sigma(k)` implies that `k = 2^{0,1} * a^2` "Squares and twice squares" [A028982](https://oeis.org/A028982).
From the sequence definitian `sigma(k + odd) == sigma(k) + odd`. So one of `sigma(k + odd)` or `sigma(k)` is odd.
So all we need to look at is terms near squares (and twice squares).
Have to look both above squares, `(16, 16+17)`, and below, `(6+19, 25)`.

* DIST <= 23, checked up to 10^18 (111 minutes)
  * Checked up to 10^17 with `pypy check_squares.py` (took 480 minutes, c++ takes 32 minutes)
* DIST <= 99, checked up to 10^17 (116 minutes)
  * See notes above for the few small terms found
* DIST <= 200000, checked up to 10^16
