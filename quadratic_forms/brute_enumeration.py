#!/usr/bin/env python

import itertools

from collections import defaultdict, Counter
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

            if temp > 0:
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
    primes = [2]
    for p in range(3, n+1, 2):
        if not any(p % q == 0 for q in primes):
            primes.append(p)
    print("primes:", print_set_short_repr(primes))

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

        common_mods = tuple(sorted(m for m, c in Counter(p % modulo for p in primes).items() if c > 3))
        extra_mods = tuple(sorted(set(mods) - set(common_mods)))
        print("\t", modulo, "\tmods", common_mods, " + ", extra_mods)

        # Exclude common mods
        test_pop = merge_primes([p for p in primes if p % modulo not in extra_mods], n)
        yield test_pop, f"exclude p % {modulo} not in {extra_mods}"

        # Exclude modulo
        if modulo in primes:
            test_pop = merge_primes([p for p in primes if p != modulo], n)
            yield test_pop, f"p != {modulo}"


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

            test_pop = merge_primes([a**2 for a in a_primes] + b_primes, n)
            yield test_pop, f"p % {modulo} in {m_split} to an even power"

            #test_pop = merge_primes([a ** 2 for a in a_primes] + [b ** 2 for b in b_primes], n)
            #yield test_pop, f"p % {modulo} all even power"



def find_brute(n = 1000):
    names = {
        (1, -2): "A035251 / A000047",
        (1, 1): "A001481 / A000050",
        (1, 2): "A002479 / A000067",
        (1, 3): "A003136 / A000205",
        (5, -1): "A031363",
    }

    min_n = 20

    forms = defaultdict(str)

    # TODO: Negative sequences are hard to enumerate
    # Find some general trick?
    for a in range(1, 100+1):
        for b in range(a, 100+1):
            population = enumerate_quadratic_form(a, b, n)
            if not population:
                continue

            population = tuple(sorted(p for p in population if p >= min_n))
            name = names.get((a, b), f"{(a, b)}")

            other_name = forms.get(population)
            #assert other_name is None, f"{name} same as {other_name}"
            if other_name:
                print(f"{name} same as {other_name}")

            forms[population] += " " + name

    print(f"Enumerated {len(forms)} Quadratic Forms")
    print()
    for test, split_name in enumerate_prime_splits(n):
        test = tuple(sorted(p for p in test if p >= min_n))
        if test in forms:
            print()
            print("Match!", split_name, "\t", forms[test])
            print("\t", print_set_short_repr(test, 24))
            print()



if __name__ == "__main__":
    find_brute(2000)
