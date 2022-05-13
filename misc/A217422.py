from collections import Counter
import itertools
import primesieve
import string

primes = primesieve.primes(81)

max_n = 10 ** 100

found = []
for r in itertools.count():
    two_s = 2 ** r
    if two_s > max_n:
        break

    n = two_s
    for s in itertools.count():
        t = str(n)
        if all(d in t for d in string.digits) and all(c in primes for c in set(Counter(t).values())):
            found.append(n)
            print ("\t", r, s, n, set(Counter(t).values()))

        n *= 17
        if n > max_n:
            break

print()
for n, an in enumerate(sorted(found)):
    print (n, an)
