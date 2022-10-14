# [A217259](https://oeis.org/A217259)

Should maybe be moved under [A050507](https://oeis.org/A050507)

This short line laid down a challange that cost me two nights of sleep.

> No term found below 2 \* 10^9 to continue sequence 435, 8576, 8826, ... - Michel Marcus, Mar 19 2013

code has now been tested up to `3e12 = 3,000,000,000,000`

Because `Sigma(i) = i + 1` iff i is prime, `is_prime(i) = (Sigma(i) - i == 1)`
Can use `sigma[i]` to check primality.


The python code can be used to search for terms in [A054903](https://oeis.org/A054903)
