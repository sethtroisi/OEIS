
# Trying to understand the work to find all
# (p, r) r = rho(p) = order2(p) = first 2^r-1 that p divides
# for p < 2^64 / r

# small r (r <= ~1200 are fully factored)
#   http://www.factordb.com/index.php?query=2%5En-1&use=n&n=1200&VP=on&VC=on&EV=on&OD=on&PRP=on&CF=on&U=on&C=on&perpage=20&format=1&sent=Show
#   (all known up to 1207 but probably much further in practice)
#   https://www.mersenne.ca/prp.php?show=0&min_exponent=1&max_exponent=2%2C000


# for large r (e.g. 12345) where it's unknown if small unknown factor exists
# have to check all p up to ~= x = 2^64 / r
# probably only have to check p = a * (factor of r) + 1
#   when r is even this is a huge number of numbers

import re

import MathLib

with open("a001265.txt") as f:
    data = f.readlines()

res = [0,0]
for row in data:
    row = row.strip()
    if row.startswith("#") or not row:
        continue

    if 'C' in row:
        print ("ERROR:", row)
        continue

    n, *factors = map(int, re.split(r'[: .]+', row))

#    print(n, f)
#    print()

    if n <= 100:
        continue

    nfactors = MathLib.factor(n)

    for f in factors:
        some = any(f % nf == 1 for nf in nfactors)
        assert some, (n, nfactors, f, factors)

