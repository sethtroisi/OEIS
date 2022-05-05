#!/usr/bin/env python

import itertools

from collections import defaultdict, Counter
from typing import List, Set, Tuple


def gen_primes(n):
    primes = [2]
    for p in range(3, n+1, 2):
        if not any(p % q == 0 for q in primes):
            primes.append(p)
    return primes

def print_set_short_repr(l, count=10) -> str:
    if len(l) <= count:
        return str(l)
    return " ".join(map(str, l[:count//2])) + " ... " + " ".join(map(str, l[-count//2:]))


def enumerate_quadratic_form(a: int, b: int, c: int, n: int) -> Set[int]:
    """Enumerate all numbers a x^2 + b y^2 <= n."""
    population = set()
    stop = n if c >= 0 else abs(c) * n

    for x in range(n):
        temp_x = a * x*x
        if temp_x > n:
            break

        for y in range(n):
            temp = temp_x + c*x*y +  b * y*y
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
        if prime == 1:
            continue
        nextGen = gen[:]
        primePower = prime
        while primePower <= n:
            for i in gen:
                number = i * primePower
                if number > n:
                    break
                nextGen.append(number)
            primePower *= prime
        gen = sorted(nextGen)
    #gen.remove(1)
    return tuple(gen)


def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s)+1))


def raise_to_power(l: List[int], p: int) -> List[int]:
    """Raise each item in l to p"""
    return tuple([a ** p for a in l])


def enumerate_prime_splits(n: int):
    primes = gen_primes(n)
    print("primes:", print_set_short_repr(primes))

    prime_splits = set()

    # Exclude 2
    test_pop = merge_primes([p for p in primes if p > 2], n)
    yield test_pop, "p > 2"

    # Exclude 2/3
    test_pop = merge_primes([p for p in primes if p > 3], n)
    yield test_pop, "p > 3"

    # TODO how to enumerate more / better splits
    for modulo in range(3, 100+1):
        mods = tuple(sorted(set(p % modulo for p in primes)))
        # TODO maybe loosen this later
        if len(mods) > 13:
            continue

        common_mods = tuple(sorted(m for m, c in Counter(p % modulo for p in primes).items() if c > 3))
        extra_mods = tuple(sorted(set(mods) - set(common_mods)))
        print("\t", modulo, "\tmods", common_mods, " + ", extra_mods)

        # Exclude common mods
        test_pop = merge_primes([p for p in primes if p % modulo not in extra_mods], n)
        yield test_pop, f"exclude p % {modulo} not in {extra_mods}"

        # Exclude a single prime (modulo)
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
            a_primes = tuple(a_primes)
            b_primes = tuple(b_primes)

            if a_primes in prime_splits:
                continue
            prime_splits.add(a_primes)

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

            test_pop = merge_primes(raise_to_power(a_primes, 2) + b_primes, n)
            yield test_pop, f"p % {modulo} in {m_split} to an even power"

            if a_primes[0] ** 3 <= n and b_primes[0] ** 2 <= n:
                test_pop = merge_primes(raise_to_power(a_primes, 3) + b_primes, n)
                yield test_pop, f"p % {modulo} in {m_split} higher powers (3, 1)"

                test_pop = merge_primes(raise_to_power(a_primes, 3) + raise_to_power(b_primes, 2), n)
                yield test_pop, f"p % {modulo} in {m_split} higher powers (3, 2)"

            if a_primes[0] ** 4 <= n:
                test_pop = merge_primes(raise_to_power(a_primes, 4) + b_primes, n)
                yield test_pop, f"p % {modulo} in {m_split} higher powers (4, 1)"

                test_pop = merge_primes(raise_to_power(a_primes, 4) + raise_to_power(b_primes, 2), n)
                yield test_pop, f"p % {modulo} in {m_split} higher powers (4, 1)"


            #test_pop = merge_primes([a ** 2 for a in a_primes] + [b ** 2 for b in b_primes], n)
            #yield test_pop, f"p % {modulo} all even power"



def find_brute(n = 1000):


    """
    Can find some of these with searches like
        "even exponent" "prime" "quadratic form"
        "even power" "prime" "quadratic form"
        "quadratic form" "prime" "representation"
    """

    names = {
        (1, -2, 0): "A035251 / A000047",
        (1, 1, 0): "A001481 / A000050",
        (1, 2, 0): "A002479 / A000067",

        (1, 3, 0): "A003136 / A000205",
        (1, 1, -1): "A003136 / A000205",
        (1, 1, 1): "A003136 / A000205",

        (5, -1, 0): "A031363",
        (1, 1, 3): "A031363",
    }

    min_n = 20

    forms = defaultdict(list)

    # TODO: Negative sequences are hard to enumerate
    # Find some general trick?
    for a in range(1, 100+1):
        for b in range(a, 100+1):
            for cross in (-1, 0, 1, 2, 3):

                # Repeated forms not needed
                if (a, b, cross) in [(1, 1, 1), (1, 2, 2), (1, 1, -1), (1, 3, 3), (2, 2, -1)]:
                    continue

                population = enumerate_quadratic_form(a, b, cross, n)
                if not population:
                    continue

                population = tuple(sorted(p for p in population if p >= min_n))
                default_name = f"{a} x^2 {cross:+} xy + {b} y^2" if cross != 0 else f"{a} x^2 + {b} y^2"
                name = names.get((a, b, cross), default_name)

                other_name = forms.get(population)
                #assert other_name is None, f"{name} same as {other_name}"
                if other_name:
                    print(f"\t{name} same as {other_name}")

                forms[population].append(name)

    print(f"Enumerated {len(forms)} Quadratic Forms")
    print()
    for test, split_name in enumerate_prime_splits(n):
        test = tuple(sorted(p for p in test if p >= min_n))
        if test in forms:
            print()
            print("Match!", split_name, "  <==>  ", " | ".join(forms[test]))
            print("\tsequence:", print_set_short_repr(test, 24))
            if not any("A" in seq for seq in forms[test]):
                print("\t", test)
            print()


def subgroups():
    comb_with_repl = itertools.combinations_with_replacement

    primes = gen_primes(1000)
    Zp = 9

    def is_closed(subgroup):
        return all(a * b % Zp in subgroup for a, b in itertools.product(subgroup, repeat=2))

    all_mods = tuple(sorted(set(p % Zp for p in primes) - {0}))
    minus_two = tuple(m for m in all_mods if m != 2)
    minus_three = tuple(m for m in all_mods if m != 3)
    minus_two_three = tuple(m for m in all_mods if m > 3)

    for mods in [all_mods, minus_two, minus_three, minus_two_three]:
        missing = sorted(set(all_mods) - set(mods))

        sub_groups = [g for g in powerset(mods) if is_closed(g)]
        # empty sub_group is kinda boring
        sub_groups.remove(tuple())

        for g1, g2 in itertools.combinations(sub_groups, 2):
            if set(g1) | set(g2) == set(mods) and (not (set(g1) & set(g2))):
                print("\tA", g1, g2)


        for g in powerset(mods):
            if not g: continue

            inv = tuple(sorted(set(mods) - set(g)))
            # Check if g^2 -> !g
            if all(a * b % Zp in inv for a, b in comb_with_repl(g, 2)):
                if inv in sub_groups:
                    g_a = tuple([p for p in primes if p % Zp in g])
                    g_b = tuple([p for p in primes if p % Zp in inv])
                    m_p = tuple([p for p in primes if p % Zp in missing])
                    print("\tA^2 -> B, B^2 -> B:", g, inv, "\tmissing:", m_p)
                    print()

                    # Could try powerset(m_p) if it contains more than one
                    print("merged a^2          :", merge_primes(raise_to_power(g_a, 2) + g_b + m_p, 100))
                    print("merged b^2          :", merge_primes(g_a + m_p + raise_to_power(g_b, 2), 100))
                    if m_p:
                        print("merged (a+missing)^2:", merge_primes(raise_to_power(g_a + m_p, 2) + g_b, 100))
                        print("merged (b+missing)^2:", merge_primes(g_a + raise_to_power(g_b + m_p, 2), 100))
                    print()
                    print()


def find_quadratic_with_primes():
    n = 5
    present = [11, 19, 31, 59, 71, 79, 131, 139, 151, 179, 191, 199]
    not_present = [29, 41, 61, 89, 101, 109, 149, 181]

    present = set(present)
    not_present = set(not_present)
    assert not (present & not_present)

    # TODO: Negative sequences are hard to enumerate
    # Find some general trick?
    for a in range(-100, 100+1):
        for b in range(-100, 100+1):
            for cross in range(-100, 100+1):
                if cross*cross - 4*a*b != -4 * n:
                    continue
                print(a, cross, b)

                population = enumerate_quadratic_form(a, b, cross, max(max(present), max(not_present)))
                population = set(population)
                if (population & present) == set(present) and len(population | not_present) == 0:
                    print("Found!", (a,b,cross))


if __name__ == "__main__":
    #find_brute(8000)
    #subgroups()
    find_quadratic_with_primes()
