import itertools
import sympy
from sympy.core.power import isqrt

def count_brute(n):
    if n > 10 ** 8: return None

    count = 0
    for y in range(isqrt(n)+2):
        y2 = y*y
        for x in range(y+1):
            if y2 + y*x + x*x == n:
                count += 1

    return count


def count_circle(n):
    # Circle walk    x^2 + y^2  +  x*y
    if n > (22 * 10 ** 7) ** 2:
        return None

    count = 0

    # 0 <= x <= y
    # Largest value  is at (0, y)  => y^2
    # Smallest value is at (y, y)  => 3*y^2
    minY = isqrt(n // 3)
    maxY = isqrt(n)

    x = 0
    for y in reversed(range(minY, maxY+1)):
        test = x*x + x*y + y*y
        while test < n and x < y:
            x += 1
            test = x*x + x*y + y*y
        if test == n:
            count += 1
#            print("\t", test, " ", x, y)

    return count


def prod(prime_powers):
    """Product of list of [(prime, power), (prime, power), ...]"""
    return sympy.prod([p**power for p, power in prime_powers])


def gen_small(signature, count, primes):
    # Generate some small [est?] numbers that match the signature.
    needed = len(signature)
    smallest = prod(zip(primes[:needed], signature))

    nums = set()

    # Allow with larger primes
    for to_use in itertools.combinations(primes[:needed+3], needed):
        nums.add(prod(zip(to_use, signature)))

    # Allow with reordering primes
    for to_use in itertools.permutations(primes[:needed+2], needed):
        nums.add(prod(zip(to_use, signature)))

    return sorted(nums)[:count]


def gen_all_signatures(N, small_primes):
    signatures = {}
    for count in range(1, len(small_primes)):
        signature = (1,) * count
        smallest = sympy.prod(small_primes[:count])
        if smallest > N:
            break
        signatures[smallest] = signature
    assert smallest > N

    while signatures:
        smallest = min(signatures)
        sig = signatures.pop(smallest)
        yield sig, smallest

        # Increment any of the powers in signature
        for power in set(sig):
            # Replace an instance of power with power+1
            index = sig.index(power)
            new_signature = sig[:index] + (power + 1,) + sig[index+1:]
            new_signature = tuple(sorted(new_signature, reverse=True))

            smallest = prod(zip(small_primes, new_signature))
            if smallest <= N:
                signatures[smallest] = new_signature


def generate_signatures(N):
    # Generate a signature and several small examples
    small_primes = list(p for p in sympy.primerange(2, 200) if p % 6 == 1)

    verify_count = 2
    A = {}

    # A19 is a high point, search just past it
    for signature, smallest in gen_all_signatures(N, small_primes):
        other_small = gen_small(signature, verify_count, small_primes)
        count_ways = list(map(count_circle, other_small))
        print (f"\t{smallest:<10} {str(signature):20}", count_ways, other_small)
#        for other, ways in zip(other_small, count_ways):
#            print ("\t\t", ways, other, sympy.factorint(other))

        ways = count_ways[0]
        # Make sure we actually computed a count
        assert ways is not None

        # Verify everything with same signature has the same counts
        for count in count_ways:
            assert count in (ways, None)

        if ways not in A:
            A[ways] = smallest
        else:
            assert A[ways] <= smallest

    # TODO what to do about missing ranges?
    last = 0
    for m, smallest in sorted(A.items()):
        if last + 2 == m:
            print(f"{last+1} ??? (> {N})")
        elif last + 2 < m:
            print(f"{last+1}-{m-1} ??? (> {N})")
        print (m, smallest)
        last = m


def probe_pattern():
    small_primes = list(p for p in sympy.primerange(2, 200) if p % 6 == 1)

    groups = [
        [(i,) for i in range(1,16)],
        [(1,) * i for i in range(1,7)],
        [(2,) * i for i in range(1,5)],
        [(3,) * i for i in range(1,5)],
        [(4,) * i for i in range(1,5)],
        [(2,) + (1,) * i for i in range(6)],
        [(2,2) + (1,) * i for i in range(5)],
        [(3,) + (1,) * i for i in range(6)],
        [(3,3) + (1,) * i for i in range(5)],
        [(3,2) + (1,) * i for i in range(5)],
        [(4,) + (1,) * i for i in range(6)],
        [(4,4) + (1,) * i for i in range(5)],
        [(4,2) + (1,) * i for i in range(5)],
        [(4,3) + (1,) * i for i in range(5)],

        [(3,), (3,3), (3,3,2), (3,3,2,2), (3,3,2,2,1)],
        [(2,), (2,2), (2,2,3), (2,2,3,3),],
    ]

    for group in groups:
        for signature in group:
            smallest = gen_small(signature, 1, small_primes)[0]
            count = count_circle(smallest)
            print (f"\t{smallest:<15} {str(signature):20} {count}")
        print()

    # (i,)                  => i//2 + 1
    # (a, b, c, 1)          => 2 * (a, b, c)
    # (a, b, c, 1, 1)       => 4 * (a, b, c)
    # (a, b, c) + n*(1,)    => 2^n * (a, b, c)
    # (a, b, c, 3) + n*(3,) => 4^n * (a, b, c, 3)
    # (a, b, c, 3) + n*(3,) => 4^n * (a, b, c, 3)

if __name__ == "__main__":
    A13 = 68574961
    #A19 = 21169376772835837
    # This will generate terms at-least through A13/A19
    generate_signatures(A19)

    # To help understand how to map signature to count
    probe_pattern()



