"""Find number of a quadratic form, Factor them, guess at prime rules"""
import primesieve
import sympy

def brute(a, b, limit=2**20):
    """Find numbers of the form a*x^2 + b*y^2"""
    nums = set()
    for x in range(100):
        for y in range(x == 0, 100):
            nums.add(a * x * x + b * y * y)

    return sorted(nums)


def build(primes_p, primes_q, limit=10000):
    """
    Build numbers of the form [p][q]^2
    """
    def expand(primes):
        t = []
        def recursive(n, i):
            if i > limit:
                return
            t.append(n)
            for pi in range(i, len(primes)):
                p = primes[pi]
                tn = n * p
                if tn > limit:
                    break
                recursive(tn, pi)

        recursive(1, 0)
        return sorted(t)

    def merge(A, B):
        m = []
        for a in A:
            for b in B:
                n = a * b
                if n > limit:
                    break
                m.append(n)
        m.sort()
        return m

    partial_p = expand(primes_p)
    partial_q = expand([p*p for p in primes_q])
    print("\tprimes_p:", primes_p[:10], "->", partial_p[:10])
    print("\tprimes_q:", primes_q[:10], "->", partial_q[:10])
    found = merge(partial_p, partial_q)
    print("\tfound:", len(found), found[:20])
    found2 = merge(found, [2 ** i for i in range(30)])
    print("\tfound2:", len(found2), found2[:20])
    return found, found2


def factor(nums):
    # Primes found to even and odd powers
    primes_even = set()
    primes_odd = set()
    for n in nums:
        if n <= 2: continue

        f = sympy.factorint(n)
        for p, e in f.items():
            # p is a weird prime
            if p == 2: continue

            if e % 2 == 0:
                primes_even.add(p)
            else:
                primes_odd.add(p)

    both = primes_even & primes_odd
    xor = primes_even ^ primes_odd

    print("Nums:\t", sorted(nums)[:25])
    print("Even:\t", sorted(primes_even)[:20])
    print("Odd:\t", sorted(primes_odd)[:20])
    print("&:\t", sorted(both)[:20])
    print("^:\t", sorted(xor)[:20])

    for mod in [4, 6, 8]:
        print(f"  Even % {mod}:\t", sorted(set(p % mod for p in primes_even)))
        print(f"  Odd % {mod}:\t", sorted(set(p % mod for p in primes_odd)))
        print("  &:\t\t", sorted(set(p % mod for p in both))[:20])
        print("  ^:\t\t", sorted(set(p % mod for p in both))[:20])


def run(a, b):
    nums = brute(a, b)
    factor(nums)

special_forms = {
    "A000047": (1, -2), # Known to have special form %8 = {1,7} | {3,5}
    "A000050": (1, 1), # Known to have special form %4 = 1 | 3
    "A000067": (1, 2), # Known to have special form %8 = {1,3} | {5,7}
    "A000205": (1, 3), # Known to have special form %6 = {1,3} | {2,5}
}

forms = {
    "A000018": (1, 16),  # x^2 + 16 y^2
    "A000021": (1, 12),
    "A000024": (1, 10),
    "A000049": (3, 4),
    "A000072": (1, 4),
    #"A000074": (1, 4),   # Odds only
    "A000075": (2, 3),
    "A000076": (4, 4, 5), # 4 x^2 + 4 x y + 5 y^2
    "A000077": (1, 6),
    "A000286": (2, 5),
    "A054150": (1, 5),
    "A054151": (1, 7),
    "A054152": (1, 8),
    "A054153": (1, 9),
    "A054157": (2, 7),
    "A054159": (2, 9),
    "A054161": (3, 3),
    "A054162": (3, 5),
    "A054163": (3, 6),
    "A054164": (3, 7),
    "A054165": (3, 8),
    "A054166": (3, 9),
    "A054167": (3, 10),
    "A054169": (4, 5),
    "A054171": (4, 7),
    "A054173": (4, 9),
    "A054175": (5, 5),
    "A054176": (5, 6),
    "A054177": (5, 7),
    "A054178": (5, 8),
    "A054179": (5, 9),
    "A054180": (5, 10),
    "A054182": (6, 7),
    "A054184": (6, 9),
    "A054186": (7, 7),
    "A054187": (7, 8),
    "A054188": (7, 9),
    "A054189": (7, 10),
    "A054191": (8, 9),
    "A054193": (9, 9),
    "A054194": (9, 10),
}

if 1:
    test = special_forms # | forms
    for seq, form in special_forms.items():
        print(seq, form)
        if len(form) == 2:
            run(*form)
        print("\n")
        break


if 0:
    for e in range(2, 10):
        L = 2 ** e
        primes = primesieve.primes(L)

        # A000047
        #p_a = [p for p in primes if p % 8 in (1, 7)]
        #p_b = [p for p in primes if p % 8 in (3, 5)]

        # A000050
        #p_a = [p for p in primes if p % 4 == 1]
        #p_b = [p for p in primes if p % 4 == 3]

        # A000067
        #p_a = [p for p in primes if p % 8 in (1, 3)]
        #p_b = [p for p in primes if p % 8 in (5, 7)]

        # A000205 result is found1 not found2
        #p_a = [p for p in primes if p % 6 in (1, 3)]
        #p_b = [p for p in primes if p % 6 in (2, 5)]


        f1, f2 = build(p_a, p_b, limit=L)
        print(e, len(f1), len(f2))

