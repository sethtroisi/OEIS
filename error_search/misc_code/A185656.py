import primesieve
import itertools

for n in range(2, 12+1):
    it = primesieve.Iterator()
    divisor = 7 ** n
    smarandache_wellin_mod = 0
    while True:
        prime = it.next_prime()
        smarandache_wellin_mod *= 10 ** len(str(prime))
        smarandache_wellin_mod += prime
        smarandache_wellin_mod %= divisor

        if smarandache_wellin_mod == 0:
            print(n, prime)
            break
