import sys
import primesieve
import gmpy2

assert len(sys.argv) == 2, "Must pass STOP"

STOP = int(sys.argv[1])
print("Running up to", STOP)

assert 11 in primesieve.primes(11)

counts = [0 for i in range(12)]
for p in primesieve.primes(STOP):
    counts[0] += 1
    for d in (2,4,6,8,10):
        if gmpy2.is_prime(p + d):
            counts[d] += 1
            if p + d > STOP:
                print("Found pair PAST stop", d, "\t", p, p + d)

print("Found", counts[0], "primes")
for d, count in enumerate(counts):
    if d == 0: continue
    if count > 1:
        print(f"\t{d} -> {count} pairs")
