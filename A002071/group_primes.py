import math
import primesieve

primes = primesieve.primes(32,200);

t = []
for p in primes:
    q = math.prod(t + [p])
    if q.bit_length() > 32:
        mult = "*".join(map(str, t))
        print(f"\t{{{mult}u, {len(t)}}},\t// {math.log2(math.prod(t)):.1f} bits")
        t = [p]
    else:
        t.append(p)
