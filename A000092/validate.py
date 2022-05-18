#!/usr/bin/env python

"""
Validate P(n) (from bfiles & README) using high precision rounding
"""

import math
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

    I, N, Pn, An = [], [], [], []
    with open("README.md") as f:
        for l in f.readlines():
            if "n = A000099" in l:
                break
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


def validate(N, An, Pn):
    getcontext().prec = 50

    M_PI = pi()
    MULT = M_PI * 4 / 3

    def V(n):
        return MULT * Decimal(n ** 3).sqrt()

    all_good = True
    for n, an, pn in zip(N, An, Pn):
        vn = V(n)
        pn_test = an - vn
        rounded = pn_test.to_integral()
        if abs(rounded) != pn:
            all_good = False
            print(f"ERROR with rounding! {n=}, {an=}")
            print(f"\t{pn} vs {pn_test}")

    return all_good

for n, an, pn in get_data():
    assert len(n) == len(an) == len(pn)
    print("Validating", len(n), "entries")
    if not validate(n, an, pn):
        print()
        print("Had Error")
