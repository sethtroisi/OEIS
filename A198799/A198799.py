import itertools
import functools

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


def count_signature(signature):
    solutions = sympy.prod((exp+1) for exp in signature)
    return (solutions + 1) // 2


def verify_pattern():
    small_primes = list(p for p in sympy.primerange(2, 200) if p % 6 == 1)

    groups = [
        [(i,) for i in range(1,16)],
        [(1,) * i for i in range(1,7)],
        [(2,) * i for i in range(1,5)],
        [(3,) * i for i in range(1,5)],
        [(4,) * i for i in range(1,4)],
        [(2,) + (1,) * i for i in range(6)],
        [(2,2) + (1,) * i for i in range(5)],
        [(2,2,2) + (1,) * i for i in range(4)],
        [(3,) + (1,) * i for i in range(6)],
        [(3,3) + (1,) * i for i in range(5)],
        [(3,3,3) + (1,) * i for i in range(4)],
        [(3,2) + (1,) * i for i in range(5)],
        [(4,) + (1,) * i for i in range(6)],
        [(4,4) + (1,) * i for i in range(5)],
        [(4,2) + (1,) * i for i in range(5)],
        [(4,3) + (1,) * i for i in range(5)],
        [(4,3,2) + (1,) * i for i in range(4)],

        [(3,), (3,3), (3,3,2), (3,3,2,2), (3,3,2,2,1)],
        [(2,), (2,2), (2,2,3), (2,2,3,3),],
    ]

    for group in groups:
        for signature in group:
            smallest = gen_small(signature, 1, small_primes)[0]
            count = count_circle(smallest)
            count2 = count_signature(signature)
            print (f"\t{smallest:<15} {str(signature):20} {count} {count2} {'  MISMATCH' * (count != count2)}")
        print()


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
    # Determine count of ways
    small_primes = list(p for p in sympy.primerange(2, 400) if p % 6 == 1)

    tested_sigs = 0
    A = {}
    verify_count = 4

    # A19 is a high point, search just past it
    for signature, smallest in gen_all_signatures(N, small_primes):
        tested_sigs += 1

        ways = count_signature(signature)
        assert ways >= 1
        if ways not in A:
            A[ways] = smallest
        else:
            assert A[ways] <= smallest

        if tested_sigs <= 1000 or tested_sigs % 10 == 0:
            print (f"\t{smallest:<10} {str(signature):20} {ways}")

            if smallest < 1e12:
                other_small = gen_small(signature, verify_count, small_primes)
                # Filter to reasonable fast verification
                other_small = [small for small in other_small if small < 1e12]
                if other_small:
                    count_ways = list(map(count_circle, other_small))
                    print (f"\t\tVerified", count_ways, other_small)

                    # Verify numbers with same signature have the same counts
                    for count in count_ways:
                        assert count in (ways, None)

    print ()
    print ("\ttested {} up to {},\n\t{} unique m's ({} to {})".format(
        tested_sigs, N, len(A), min(A), max(A)))
    print ()

    for prime in sympy.primerange(2, 100):
        print (f"A({prime}) prod((e_i + 1)^count) in ({2*prime},{2*prime-1})")
        if prime in A:
            num = A[prime]
            print ("\t", prime, "\t", num, sympy.factorint(num))


    print ()
    last = 0
    for m, smallest in sorted(A.items()):
        if last + 2 == m:
            print(f"{last+1} ??? (> {N})")
        elif last + 2 < m:
            print(f"{last+1}-{m-1} ??? (> {N})")
        print (m, smallest)
        last = m

        # Sparse after this.
        if m > 100: break


@functools.lru_cache(maxsize=None)
def all_factorizations(n, max_d=None):
    # i.e 12 => [[12], [6*2], [4*3], [3*2*2]]

    if max_d is None:
        max_d = n

    f = []
    if n <= max_d:
        f.append((n,))

    for d in reversed(range(2, min(n, max_d+1))):
        if n % d == 0:
            for fs in all_factorizations(n // d, d):
                f.append((d,) + fs)
    return f


def factordb_format(number):
    if number < 1e10:
        return str(number)
    strN = str(number)
    length = len(strN)
    if number < 1e24:
      return "{}<{}>".format(strN, length)
    return "{}...{}<{}>".format(strN[:10], strN[-2:], length)


def smallest_m_ways(m):
    if m == 0:
        # First number not on a circle
        return 2
    if m == 1:
        # First number with only one way
        return 0

    # For the prime factorization of n
    # let S_1 be the set of distinct prime factors p_i for which p_i == 1 (mod 3),
    # let S_2 be the set of distinct prime factors p_j for which p_j == 2 (mod 3),
    # let M be the exponent of 3.
    # n = 3^M * (Prod_{p_i in S_i} p_i ^ e_i) * (Prod_{p_j in S_j} p_j ^ e_j)
    # Number of representations of n as x^2+xy+y^2=n, 0 <= x <= y is
    #   m = (Product_{p_i in S_1} (e_i + 1) + 1) // 2

    # Solving for m
    #   (prod( (e_i + 1) ) + 1) // 2
    small_primes = list(p for p in sympy.primerange(2, 200) if p % 6 == 1)

    smallest = None

    for mult in [2 * m - 1, 2 * m]:
        # factorization of mult to y pieces
        for factorization in all_factorizations(mult):
            assert len(small_primes) >= len(factorization)

            signature = [f - 1 for f in factorization]
            smallest_sig = prod(zip(small_primes, signature))
            if smallest is None or smallest_sig <= smallest:
                smallest = smallest_sig

            print ("\t{:3} {:3} {:20} {}".format(
                m, mult, str(factorization), factordb_format(smallest_sig)))

    return smallest


def gen_sequence(n):
    with open("b198799.txt", "w") as f:
        for m in range(0, n+1):
            an = smallest_m_ways(m)
            print (m, an)
            assert len(str(an)) <= 990
            f.write("{} {}\n".format(m, an))


if __name__ == "__main__":
    # Visual explination of how to map signature to count
    #verify_pattern()

    #A13 = 68574961
    #A19 = 21169376772835837

    # This will generate terms at-least through A13/A19
    #generate_signatures(A19)

    # Generate the sequence more intelligently
    gen_sequence(1000)

