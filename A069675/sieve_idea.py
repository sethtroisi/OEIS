from math import sqrt
from tqdm import tqdm

def sieveOfErat(Stop):
	if Stop < 2: return []
	length = ((Stop//2)-1 + Stop % 2)
	sieve = [True] * (length)
	for i in range(int(sqrt(Stop)) >> 1):
		if not sieve[i]: continue
		for j in range( (i*(i + 3) << 1) + 3, length, (i << 1) + 3):
			sieve[j] = False
	return [2] + [(i << 1) + 3 for i,test in enumerate(sieve) if test]

MAX = 10 ** 12
#primes = set(sieveOfErat(MAX))
primes = {a * 10 ** d + b: (d, a, b)
    for d in range(1, 30) for a in range(10) for b in [1,3,7,9]
        if (a,b) != (1,1) and (a + b) % 3 != 0}


combos = (
    ([1, 6], [24, {1, 7}]),
    ([1, 10], [40, {1, 9, 11, 19}]),
    ([1, 30], [120, {1, 31, 49, 79}]),
)

hits = 0
for x in tqdm(range(1, MAX)):
    x2 = x*x
    if x2 > MAX: break

    for y in range(1, MAX):
        y2 = y*y
        if y2 > MAX: break
        if x > 10 and y > 10: break


        for combo in combos:
            (xm, ym), (m, modulo) = combo
            t = xm * x2 + ym * y2
            if t > MAX: continue

            if t in primes:
                assert t % m in modulo, (x, y, combo)
                hits += 1
                if t > 100 and set(str(t)[1:-1]) == {"0"}:
                    print (t, x, y, "\t", combo, "\t", primes[t])

print(hits)
