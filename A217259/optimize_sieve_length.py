#!/usr/bin/python

"""
Search for good sieve_length (as calculated by the number of factors avoided)
"""

from collections import Counter
import math


def factor_integer(n):
    factors = []
    d = 2
    while d*d <= n:
        while n % d == 0:
            factors.append(d)
            n //= d
        d += 1
    if n > 1:
        factors.append(n)
    return factors


def gen_divisors(sieve_length):
    factors = factor_integer(sieve_length)
    primes = sorted(Counter(factors).items())

    def generate(i):
        prime, exp = primes[i]
        if i == 0:
            pp = 1
            for exp in range(exp+1):
                yield pp
                pp *= prime
            return

        powers = [prime ** exp for exp in range(exp+1)]
        for partial in generate(i-1):
            for pp in powers:
                yield partial * pp

    return sorted([d for d in generate(len(primes)-1)])


def calculate_covered(sieve_length, max_precalculated):
    """ Inverse sum of all divisors of sieve_length <= max_precalculated"""
    covered = 0
    divisors = [d for d in gen_divisors(sl) if d in range(2, max_precalculated+1)]
    for divisor in divisors:
        if divisor <= max_precalculated:
            covered += sieve_length // divisor
    return covered, divisors


assert list(gen_divisors(2*3*2*5)) == [1, 2, 3, 4, 5, 6, 10, 12, 15, 20, 30, 60]

START = 1008
MULT = 2*2*3
assert START % MULT == 0  # OTHERWISE sl can't divide certain numbers

best_covered = 0
best = 0
last_at = 0
for sl in range(1000, 1500000, 2*2):
    covered, divisors = calculate_covered(sl, 600)
    num_per = covered / sl
    if num_per > best:
        best = num_per
        last_at = sl
        print(f"\t{sl:6} covers {covered:6} = {num_per:2f} BEST", divisors[:10])
    if covered > best_covered and sl > 1.07 * last_at:
        best_covered = covered
        last_at = sl
        print(f"\t{sl:6} covers {covered:6} = {num_per:2f}")
