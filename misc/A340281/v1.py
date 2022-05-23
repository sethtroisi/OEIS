import math

import primesieve
import sympy
from tqdm import tqdm

def T(n, k):
    if k % 2 == 0:
        return 1

    count_a, count_b = 0, 0
    for m in range(n):
        z = pow(m, k, n)
        if pow(m, k, n) == m:
            count_a += 1
        if (-z % n) == m:
            count_b += 1

    if k % 2 == 0:
        assert count_a == count_b, (count_a, count_b)
    assert count_a % count_b == 0
    return count_a // count_b

def A340281_brute(n):
    it = primesieve.Iterator()
    def prime_iterator():
        while True:
            yield it.next_prime()

    for prime in prime_iterator():
        unique = len(set(T(prime, k) for k in range(prime)))
        if unique == n:
            return prime

def remove_twos(n):
    while n > 1 and (n & 1) == 0:
        n >>= 1
    return n

def A340281_jinyuan_wang(n):
    """Least prime p such that p-1 has n-1 odd divisors."""
    it = primesieve.Iterator()
    def prime_iterator():
        while True:
            yield it.next_prime()

    # Odd divisors remove all factors of 2
    for prime in prime_iterator():
        k = remove_twos(prime - 1)
        assert k % 2 != 0, (prime, k)

        factors = sympy.factorint(k)
        # check if n-1 odd divisors
        # print(n, prime, factors, math.prod(factors.values()))

        if math.prod([e + 1 for e in factors.values()]) == n-1:
            return prime

for n in range(2, 10):
    print(n, A340281_brute(n), A340281_jinyuan_wang(n))
