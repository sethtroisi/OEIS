import sys
import primesieve
import gmpy2
import sympy

assert len(sys.argv) == 2, "Must pass STOP"

STOP = int(sys.argv[1]) - 1
MAX_DIST = 12

print("Running up to", STOP, "checking up to", MAX_DIST-1, "away")
print()


print("Checking near STOP")
sigmas = {a: sympy.ntheory.factor_.divisor_sigma(a) for a in range(STOP-2, STOP+MAX_DIST+1)}
for a in range(min(sigmas), STOP+1):
    sigma_a = sigmas[a]
    for d in range(2, MAX_DIST):
        if (a + d) in sigmas and (sigma_a + d == sigmas[a + d]):
            print(f"Match near STOP, sigma({a}) = sigma({a + d}) = {sigma_a}")
print("Finished STOP checks")
print()

# Verify my understanding of primesieve endpoint inclusion
assert 11 in primesieve.primes(11)

counts = [0 for i in range(MAX_DIST+1)]
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
