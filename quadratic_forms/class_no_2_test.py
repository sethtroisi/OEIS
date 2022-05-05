import primesieve
import gmpy2
import brute_enumeration

M = 100

n = 5
print(f"{n=}")
print()
print(f"u^2 + {n} v^2")

if n in (5, 13, 37):
    a, b, c = 2, 2, (n+1)//2
    # There are many forms here with different signs, but they all generate the same numbers
    print(f"{a} u^2 + {b} uv + {c} v^2")
    print()
elif n in (6,10,22,58):
    a, b, c = 2, 0, n // 2
    print(f"{a} u^2 + {b} uv + {c} v^2")
    print()
elif n in (-83,-67,-59,-43,-19,-11,-3):
    a, b, c = -n, 0, -1
    print(f"{a} u^2 - v^2")
    print()
else:
    assert(f"Unsupported {n=}")

assert b*b - 4*a*c == -4 * n

P = []
Q = []
for p in primesieve.primes(3, M):
    l = gmpy2.legendre(-n, p)
    if l == 1:
        P.append(p)
    if l == -1:
        Q.append(p)
print(f"{P=}")
print(f"{Q=}")
print()

P0_n = sorted(set(z for z in (u*u + n * v*v for u in range(-1000, 1000) for v in range(-1000, 1000)) if 0 <= z <= M))
P1_n = sorted(set(z for z in (a*u*u + b*u*v + c*v*v for u in range(-1000,1000) for v in range(-1000,1000)) if 0 <= z <= M))
print(f"{P0_n=}")
print(f"{P1_n=}")
print()

print("for copy and pasting")
# Primes from P0_n and P1_n
P0 = [p for p in P if p in P0_n]
P1 = [p for p in P if p in P1_n]

print(f"{P0=}")
print(f"{P1=}")
print()
assert not (set(Q) & set(P0_n))


missing_P = sorted(set(P) - set(P0_n) - set(P1_n))
if missing_P:
    print("Missing P:", missing_P)
    assert not missing_P

# Need to know where 2 and p appear
prime_divisors = [2, n if n % 2 == 1 else n // 2]
assert gmpy2.is_prime(prime_divisors[1])

P0_d = [p for p in prime_divisors if p in P0_n]
P1_d = [p for p in prime_divisors if p in P1_n]
print(f"{P0_d=}")
print(f"{P1_d=}")
print()
assert len(P0_d + P1_d) == 2, P0_d

mp = brute_enumeration.merge_primes
left  = mp(P0 + P0_d + [q*q for q in Q if q*q < M], M)

# TODO how to generate even multiplicity?
right = mp(P1_d + P1, M)

all_members = sorted(set(mp(left + right, M)))
print(all_members)
