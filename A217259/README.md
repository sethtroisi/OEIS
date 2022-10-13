# [A217259](https://oeis.org/A217259)

This short line laid down a challange that cost me a night of sleep

> No term found below 2 \* 10^9 to continue sequence 435, 8576, 8826, ... - Michel Marcus, Mar 19 2013

code has now been tested up to `200,000,000,000`

Because `Sigma(i) = i + 1` iff i is prime, `is_prime(i) = (Sigma(i) - i == 1)`V
Can use `sigma[i]` to check primality.

Speed up GMP `mpz_probab_prime_p` using `is_prime` from
["Fast Primality Testing for Integers That Fit into a Machine Word"](http://ceur-ws.org/Vol-1326/020-Forisek.pdf)
Which is significantly faster than GMP (because of machine word sized inputs). But this is not needed anymore.



