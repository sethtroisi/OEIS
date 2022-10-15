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


