import array
import itertools
import math
from decimal import Decimal, getcontext
from tqdm import tqdm
import primesieve

# From https://docs.python.org/3/library/decimal.html
def pi():
    """Compute Pi to the current precision.

    >>> print(pi())
    3.141592653589793238462643383

    """
    getcontext().prec += 2  # extra digits for intermediate steps
    three = Decimal(3)      # substitute "three=3.0" for regular floats
    lasts, t, s, n, na, d, da = 0, three, 3, 1, 0, 0, 24
    while s != lasts:
        lasts = s
        n, na = n+na, na+8
        d, da = d+da, da+32
        t = (t * n) / d
        s += t
    getcontext().prec -= 2
    return +s               # unary plus applies the new precision


def get_n2_counts(N_2):
    # I wish array had something like numpy.zeros

    counts = array.array('H', [0]) * (N_2+1)
    assert counts.itemsize == 2, counts.itemsize
    assert len(counts) == N_2 + 1

    tuples = 0
    for i in tqdm(itertools.count(), total=math.isqrt(N_2) + 1):
        i_2 = i*i
        if i_2 > N_2:
            break

        # Handle i == j
        if 2 * i_2 <= N_2:
            tuples += 1
            counts[2 * i_2] += (1 + 3 * (i > 0))  # sign

        # Handle j == 0
        if i != 0:
            tuples += 1
            counts[i_2] += 2 * 2  # order * sign

        for j in range(1, i):
            i_j = i_2 + j*j
            if i_j > N_2:
                break

            tuples += 1
            counts[i_j] += 8

    print("\ttuples:", tuples)
    print("\tsum(counts):", sum(counts))
    print("\tmax(counts):", max(counts))

    return counts


def get_n2_counts_factoring(N_2):
    counts = array.array('H', [4]) * (N_2+1)
    assert counts.itemsize == 2, counts.itemsize
    assert len(counts) == N_2 + 1

    updates = 0
    counts[0] = 1

    # https://mathworld.wolfram.com/SumofSquaresFunction.html
    for p in primesieve.primes(2, N_2):
        if p == 2:
            # Special handling
            continue

        if p % 4 == 3:
            for exp in range(1, 100, 2):
                # zero any multiple not to an even power
                pp = p ** exp
                if pp > N_2:
                    break

                ppp = p * pp

                # Zero any odd power
                # any m * pp, m % p != 0
                for i_pp in range(0, N_2, ppp):
                    stop = min(i_pp + ppp, N_2 + 1)
                    for m in range(i_pp + pp, stop, pp):
                        counts[m] = 0
                        updates += 1


        if p % 4 == 1:
            # Need to handle each exponent
            for exp in range(1, 100):
                pp = p ** exp
                if pp > N_2:
                    break

                ppp = p * pp

                # any m * pp, m % p != 0
                for i_pp in range(0, N_2, ppp):
                    stop = min(i_pp + ppp, N_2 + 1)
                    for m in range(i_pp + pp, stop, pp):
                        counts[m] *= 1 + exp
                        updates += 1

    print("\tupdates:", updates)
    return counts


def get_n2_counts_on_demand_factoring(N_2):
    factor = array.array('I', [4]) * (N_2+1)
    assert factor.itemsize == 4, factor.itemsize
    assert len(factor) == N_2 + 1

    updates = 0

    r = math.isqrt(N_2)

    # https://mathworld.wolfram.com/SumofSquaresFunction.html
    for p in primesieve.primes(2, r):
        if p % 4 == 1:
            for m in range(p * p, N_2+1, p):
                factor[m] = p
                updates = 0

    for p in primesieve.primes(2, N_2):
        if p % 4 == 3:
            # time using p vs p^2 here
            for m in range(p * p, N_2+1, p):
                factor[m] = 0
                updates += 1


    counts = array.array('H', [4]) * (N_2+1)
    assert counts.itemsize == 2, counts.itemsize
    assert len(counts) == N_2 + 1
    counts[0] = 1

    for i, f in enumerate(factor[2:], 2):
        m = 4
        if f == 0:  # has factor of p
            continue

    print("\tupdates:", updates)
    return counts


N = 2 * 10 ** 7
counts1 = get_n2_counts(N)
counts2 = get_n2_counts_on_demand_factoring(N)

#for i, (c1, c2) in enumerate(zip(counts1, counts2)):
#    if c1 != c2:
#        print ("mismatch at", i, "\t", c1, c2)

assert(counts1 == counts2)


# TODO would be faster to construct signature with
# some_factor[n]
# with some_factor[n] prioritizing p % 4 == 3 factors
# then build 4 * (1 + exp) * (1 + exp)

