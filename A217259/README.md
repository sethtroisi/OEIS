# [A217259](https://oeis.org/A217259)

This short line laid down a challange that cost me a night of sleep

> No term found below 2*10^9 to continue sequence 435, 8576, 8826, ... - Michel Marcus, Mar 19 2013

code has now been tested up to `200,000,000,000`

Tried to speed up gmp `mpz_probab_prime_p` using [arduino Prime64](https://github.com/going-digital/Prime64). Required after changing `mulMod` to make use of `__uint128_t` it was 10% faster, but not worth the uncertainty.

If `Sigma(i) = i + 1` iff i is prime, I can speed up the code 2x.

