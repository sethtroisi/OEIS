import primesieve
import itertools
from typing import List, Set, Tuple


def print_set_short_repr(l, count=10) -> str:
    if len(l) <= count:
        return str(l)
    return " ".join(map(str, l[:count//2])) + " ... " + " ".join(map(str, l[-count//2:]))


def enumerate_quadratic_form(a: int, b: int, n: int) -> Set[int]:
    """Enumerate all numbers a x^2 + b y^2 <= n."""
    population = set()
    for x in range(n):
        temp_x = a * x * x
        if temp_x > n:
            break

        for y in range(n):
            temp = temp_x + b * y * y
            if temp > n:
                break
            population.add(temp)

    return population


def merge_primes(primes: List[int], n: int) -> Tuple[int, ...]:
    """
    Merge "primes" (some are prime^2)

    Similar to Smooth Numbers
    Returns all numbers ( <= n) composed only of primes from primes

    Can be speed up if so desired
    """

    gen = [1]
    for prime in primes:
        nextGen = gen[:]
        primePower = 1
        while primePower <= n:
            primePower *= prime
            for i in gen:
                number = i * primePower
                if number > n:
                    break
                nextGen.append(number)
        gen = sorted(nextGen)
    #gen.remove(1)
    return tuple(gen)


def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s)+1))


def enumerate_prime_splits(n: int):
    primes = primesieve.primes(0, n)

    prime_splits = set()

    # TODO how to enumerate more / better splits
    for modulo in range(3, 100+1):
        # TODO maybe loosen this later
        #if modulo in primes:
        #    continue

        mods = tuple(sorted(set(p % modulo for p in primes)))
        # TODO maybe loosen this later
        if len(mods) > 12:
            continue

        print("\t", modulo, "\t", mods)

        for m_split in powerset(mods):
            # Don't allow the full / empty set
            if not m_split or m_split == mods:
                continue

            # Consider itertool partition recipe
            a_primes = []
            b_primes = []
            for p in primes:
                if p % modulo in m_split:
                    a_primes.append(p)
                else:
                    b_primes.append(p)

            if tuple(a_primes) in prime_splits:
                continue
            prime_splits.add(tuple(a_primes))

            """
            Try only a primes
            Try a^1, b^2
            Inverses are tested with the inverse m_split
            """

#            print()
#            print(f"\t{m_split!s:12}", print_set_short_repr(a_primes))
#            print(f"\t{' ':12}", print_set_short_repr(b_primes))

            test_pop = merge_primes(a_primes, n)
            yield test_pop, f"p % {modulo} in {m_split}"

            test_pop = merge_primes(a_primes + [b ** 2 for b in b_primes], n)
            yield test_pop, f"p % {modulo} not in {m_split} even power"

            #test_pop = merge_primes([a ** 2 for a in a_primes] + b_primes, n)
            #yield test_pop, f"p % {modulo} in {m_split} even power"

            test_pop = merge_primes([a ** 2 for a in a_primes] + [b ** 2 for b in b_primes], n)
            yield test_pop, f"p % {modulo} all even power"



def find_brute(n = 2000):
    min_n = 20

    forms = {}
    for a in range(1, 100+1):
        for b in range(a, 100+1):
            population = enumerate_quadratic_form(a, b, n)
            if not population:
                continue

            population = tuple(sorted(p for p in population if p >= min_n))
            assert population not in forms, (a, b)
            forms[population] = (a, b)

    print(f"Enumerated {len(forms)} Quadratic Forms")
    print()
    for test, split_name in enumerate_prime_splits(n):
        test = tuple(sorted(p for p in test if p >= min_n))
        if test in forms:
            print()
            print("Match!", split_name, "\t", forms[test])
            print("\t", print_set_short_repr(test, 24))
            print()

find_brute()
