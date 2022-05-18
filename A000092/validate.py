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


def validate(N, An, Pn):
    getcontext().prec = 50

    assert len(N) > 0, (len(N), len(An), len(Pn))
    M_PI = pi()

    # Check if A000099 (An ~ 4/3*Pi*N^(3/2)) or A000092 (An ~ Pi*N)
    if abs(math.pi - An[-1] / N[-1]) < 0.2:
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

    return num_errors

for n, an, pn in get_data():
    assert len(n) == len(an) == len(pn)
    print("Validating", len(n), "entries")
    num_errors = validate(n, an, pn)
    if num_errors:
        print("Had", num_errors, "Error")
    print()
