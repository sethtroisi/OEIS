#!/usr/bin/env python

"""
Validate P(n) (from bfiles & README) using high precision rounding
"""

import math
import tqdm

from decimal import Decimal, getcontext, localcontext


# From https://docs.python.org/3/library/decimal.html
def pi():
    """Compute Pi to the current precision.

    >>> print(pi())
    3.141592653589793238462643383

    """
    with localcontext() as ctx:
        ctx.prec += 2  # extra digits for intermediate steps
        three = Decimal(3)      # substitute "three=3.0" for regular floats
        lasts, t, s, n, na, d, da = 0, three, 3, 1, 0, 0, 24
        while s != lasts:
            lasts = s
            n, na = n+na, na+8
            d, da = d+da, da+32
            t = (t * n) / d
            s += t
    return +s               # unary plus applies the new precision


def get_data():
    # Read b files and README.md
    N, Pn, An = [], [], []
    for fn, temp in [("b000092.txt", N), ("b000223.txt", Pn), ("b000413.txt", An)]:
        with open(fn) as f:
            for l in f.readlines():
                temp.append(int(l.split(" ")[1]))

        print(fn, "->", temp[:10])

    yield (N, An, Pn)

    N, Pn, An = [], [], []
    for fn, temp in [("b000099.txt", N), ("b000036.txt", Pn), ("b000323.txt", An)]:
        with open(fn) as f:
            for l in f.readlines():
                temp.append(int(l.split(" ")[1]))

        print(fn, "->", temp[:10])

    yield (N, An, Pn)

    with open("README.md") as f:
        lines = f.readlines()

    # Split into two sets of lines (before and after "### Results for A000099")
    split = min(i for i, line in enumerate(lines) if "### Results for A000099" in line)
    lines_A000092 = lines[:split]
    lines_A000099 = lines[split:]

    for group in [lines_A000092, lines_A000099]:
        I, N, Pn, An = [], [], [], []
        for l in group:
            if l.startswith("|"):
                values = list(r.strip() for r in l.split("|")[1:-1])
                if all(v.replace("-", "").isnumeric() for v in values):
                    i, n, p, a = map(int, values)
                    #print("\t", i, n, p, a)
                    I.append(i)
                    N.append(n)
                    Pn.append(p)
                    An.append(a)

        print("README.md ->", I[:10], I[-10:])
        yield (N, An, Pn)


def sign_order_mult(i, j, k):
    same = (i == j == k) + (i == j) + (j == k)
    mult = (6, 3, None, 1)[same]

    zeros = (i == 0) + (j == 0) + (k == 0)
    sign = (8,4,2,1)[zeros]
    return mult * sign


def validate_r3_An(n):
    """count i,j,k any order, any sign such that i^2 + j^2 + k^2 <= n"""
    count = 0

    wrapper = tqdm.tqdm if n > 1000 else lambda l: l
    for i in wrapper(range(math.isqrt(n)+1)):
        i_2 = i*i
        for j in range(min(i, math.isqrt(n - i_2)) + 1):
            max_k = math.isqrt(n - i_2 - j*j)

            # k = 0
            count += sign_order_mult(i, j, 0)

            # 1 <= k < j
            if j > 1 and max_k:
                count += min(max_k, j-1) * sign_order_mult(i, j, 1)

            # k == j
            if j >= 1 and max_k >= j:
                count += sign_order_mult(i, j, j)

    return count


# Validate our logic with a quick self-check
selfcheck = validate_r3_An(1000)
assert selfcheck == 132451, selfcheck


def validate(N, An, Pn):
    getcontext().prec = 50

    assert len(N) > 0, (len(N), len(An), len(Pn))
    M_PI = pi()

    # Check if A000099 (An ~ 4/3*Pi*N^(3/2)) or A000092 (An ~ Pi*N)
    isA000099 = abs(math.pi - An[-1] / N[-1]) < 0.2
    if isA000099:
        def V(n):
            return M_PI * n
    else:
        MULT = M_PI * 4 / 3
        def V(n):
            return MULT * Decimal(n ** 3).sqrt()

    num_errors = 0
    for n, an, pn in zip(N, An, Pn):
        vn = V(n)
        pn_test = an - vn
        rounded = pn_test.to_integral()
        if rounded != pn:
            num_errors += 1
            if num_errors < 10:
                print(f"ERROR with rounding! {n=}, {an=}")
                print(f"\t{pn} vs {rounded} from {pn_test}")

    if not isA000099:
        test = validate_r3_An(n)
        mismatch = an != test
        print(f"\tAn({n}) = {test} {'DOES NOT MATCH' if mismatch else 'matches'} {an}")
        num_errors += mismatch

    return num_errors


for n, an, pn in get_data():
    assert len(n) == len(an) == len(pn)
    print("Validating", len(n), "entries")
    num_errors = validate(n, an, pn)
    if num_errors:
        print("Had", num_errors, "Error")
    print()
