#!/usr/bin/env python

import math
import sympy

from collections import defaultdict

max_x = 500
k = {3*x*x + 4*y*y: (x, y) for x in range(max_x+1) for y in range(max_x+1)}

pairs = defaultdict(list)


seen = set()

for a, (x, y) in sorted(k.items()):
    if a == 0 or a > 3*max_x*max_x:
        continue

    # Remove scaled items (x*g, y*g)
    g = math.gcd(x, y)
    if g > 1:
        assert a // (g ** 2) in k, (a, g)
        continue

    f = sympy.factorint(a)

    # From (0, y)
    if a % 4 == 0 and all(e % 2 == 0 for e in f.values()):
        ty = sympy.prod(p ** (e//2) for p, e in f.items()) // 2
        assert (3 * 0**2 + 4 * ty**2) == a
        continue

    # From (x, 0)
    if a % 3 == 0 and all(e % 2 == 0 for p, e in f.items() if p != 3):
        tx = sympy.prod(p ** (e//2) for p, e in f.items())
        assert (3 * tx**2 + 4 * 0**2) == a
        continue

    # Check if we can scale an item
    if any(e % 2 == 0 for p, e in f.items() if p != 2):
        g = sympy.prod(p ** (e//2) for p, e in f.items() if p != 2 and e % 2 == 0)
        assert g > 1

        assert a // (g ** 2) in k
        tx, ty = k[a // (g ** 2)]
        assert (3 * (tx * g)**2 + 4 * (ty * g) ** 2) == a
        continue

    # Simple statement about the sequence
    if len(f) == 1:
        assert a % 12 == 7
        # prime of the form 12*k + 7
        continue

    # From (2*x, y) -> 4 * (3*x^2 + y^2) See A000050 representations
    if a % 4 == 0 and all(e % 2 == 0 for p, e in f.items() if p % 3 == 2):
        continue

    # Investigating even powers
    assert all(e % 2 == 1 for e in f.values()), (a, f)

    # Cuban primes (x^2 + xy + y^2 <=> x^2 + 3y^2) or p % 6 == 1
    assert all(p == 3 or p % 6 == 1 for p in f.keys()), (a, f)

    sig = tuple(sorted(f.items()))
    assert sig not in seen
    seen.add(sig)

    if len(f) == 2 and sum(f.values()) > 2:
        if 13 in f:
            a, b = f.keys()
            pairs[a].append(b)
            pairs[b].append(a)
            print(f)

    if len(f) == 2 and f.get(3) == 1:
        assert max(f) % 12 == 1
        continue

#    if len(f) == 2 and f.get(13) == 1:
#        print(f"{a:6}  {len(seen)}th from {f!s:24}   {(x, y)}")



# Any prime can occur to an even power
# (x, y) -> z | (x*p, y*p) -> z*p^2


# Any prime z must be of the form 12*k + 7
#print("odd primes:", sorted(odd_primes)[:30], "...")
#assert all (o == 3 or o % 6 == 1 for o in odd_primes)

# For two primes both to 1 power
# {3,7,19,31,43}  +  {p = 1 mod 12}
# {13,37,61,73} + ({3} + (p = 7 mod 12))

for p in sorted(pairs):
    if p < 200 and len(pairs[p]) > 1:
        print(p, p %12, pairs[p][1] % 12, "\t", sorted(pairs[p][:20]))
        #others = [o % 12 for o in pairs[p][2:]]
        #assert len(set(others)) == 1, others


