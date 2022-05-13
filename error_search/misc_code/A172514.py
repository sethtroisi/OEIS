import array
import gmpy2
import math
import primesieve
import time

"""
with b = base-1

looking for p = d1,d2,d3,d4,..dn

where none of these are prime
    1 d1 d2 ... dn 0
    1 d1 d2 ... dn 1
    ...
    1 d1 d2 ... dn b
    2 d1 d2 ... dn 0
    2 d1 d2 ... dn 1
    ...
    2 d1 d2 ... dn b
    ...
    b d1 b2 ... dn 0
    b d1 b2 ... dn 1
    ...
    b d1 b2 ... dn b

Have b different rows of sieving starting at
    i * base ** (n+2) + p * base + 0

Hard value is base = 26
L1 is to small, so size each run for L2 / base

Unfortanetly most of each sieve is unused
    (p * base) skips huge chunks of the array

Could do b passes (for each lead) over say 1000 primes at a time
    sieve[pi][trail]

For all pi with no primes (try again with lead=2)
    sieve[pi_filtered][trail]

[trail] is kinda silly because only 6-8 values are possibly prime
    wheel30 eliminates 22/30 values

...

"""

for base in range(2, 26+1):
    T0 = time.time()

    digits = 0
    next_digit = 0

    it = primesieve.Iterator()
    while True:
        prime = it.next_prime()
        if prime > next_digit:
            digits += 1
            place = base ** (digits + 1)
            #pairs = [(l, t) for l in range(1, base) for t in range(base)
            #    if math.gcd(l * place + t, base) == 1 and (prime == 2 or (l * place + base + t) % 2 == 1)]
            #print("\t", base, digits, pairs)
            next_digit = base ** digits

        num = prime * base
        found = False
        for lead_digit in range(1, base):
            temp = lead_digit * place + num
            for trail_digit in range(1 - (temp % 2), base, 2):
                n = temp + trail_digit
        #for l, t in pairs:
        #        n = l * place + num + t
                if gmpy2.is_prime(n):
                    found = True
                    break
            if found:
                break

        if not found:
            print (base, prime, f"\t{time.time() - T0:.1f}")
            break
