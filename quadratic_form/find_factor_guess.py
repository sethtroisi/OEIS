"""Find number of a quadratic form, Factor them, guess at prime rules"""
import sympy

def brute(a, b, limit=2**20):
    """Find numbers of the form a*x^2 + b*y^2"""
    nums = set()
    for x in range(100):
        for y in range(x == 0, 100):
            nums.add(a * x * x + b * y * y)

    return sorted(nums)


def factor(nums):
    # Primes found to even and odd powers
    primes_even = set()
    primes_odd = set()
    for n in nums:
        f = sympy.factorint(n)
        for p, e in f.items():
            if e % 2 == 0:
                primes_even.add(p)
            else:
                primes_odd.add(p)

    both = primes_even & primes_odd
    xor = primes_even ^ primes_odd

    print("Even:\t", sorted(primes_even)[:20])
    print("Odd:\t", sorted(primes_odd)[:20])
    print("&:\t", sorted(both)[:20])
    print("^:\t", sorted(xor)[:20])

    for mod in [6, 8]:
        print(f"  Even % {mod}:\t", sorted(set(p % mod for p in primes_even)))
        print(f"  Odd % {mod}:\t", sorted(set(p % mod for p in primes_odd)))
        print("  &:\t\t", sorted(set(p % mod for p in both))[:20])
        print("  ^:\t\t", sorted(set(p % mod for p in both))[:20])


def run(a, b):
    nums = brute(a, b)
    factor(nums)

forms = {
    "A000018": (1, 16),
    "A000021": (1, 12),
    "A000024": (1, 10),  # x^2 + 10 y^2
    "A000047": (1, -2),  # x^2 - 2 y^2
    "A000049": (3, 4),   # 3 x^2 + 4 y^2
    "A000050": (1, 1),
    "A000067": (1, 2),
    "A000072": (1, 4),
    #"A000074": (1, 4),   # Odds only
    "A000075": (2, 3),
    "A000076": (4, 4, 5), # 4 x^2 + 4 x y + 5 y^2
    "A000077": (1, 6),
    "A000205": (1, 3),
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

for seq, form in forms.items():
    print(seq, form)
    if len(form) == 2:
        run(*form)
    print("\n")
