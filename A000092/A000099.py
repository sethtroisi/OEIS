import array
import itertools
import math
import time
from collections import Counter
from decimal import Decimal, getcontext, localcontext
from tqdm import tqdm

# From https://docs.python.org/3/library/decimal.html
def pi():
    """Compute Pi to the current precision.

    >>> print(pi())
    3.141592653589793238462643383

    """
    with localcontext() as ctx:
      ctx.prec += 2         # extra digits for intermediate steps
      three = Decimal(3)    # substitute "three=3.0" for regular floats
      lasts, t, s, n, na, d, da = 0, three, 3, 1, 0, 0, 24
      while s != lasts:
          lasts = s
          n, na = n+na, na+8
          d, da = d+da, da+32
          t = (t * n) / d
          s += t
    return +s               # unary plus applies the new precision


def get_n2_counts(N_2):
    # I wish array had something like numpy.zeros

    counts = array.array('H', [0]) * (N_2+1)
    assert counts.itemsize == 2, counts.itemsize
    assert len(counts) == N_2 + 1

    tuples = 0
    for i in tqdm(itertools.count(), total=math.isqrt(N_2) + 1):
        i_2 = i*i
        if i_2 > N_2:
            break

        # Handle i == j
        if 2 * i_2 <= N_2:
            tuples += 1
            counts[2 * i_2] += (1 + 3 * (i > 0))  # sign

        # Handle j == 0
        if i != 0:
            tuples += 1
            counts[i_2] += 2 * 2  # order * sign

        for j in range(1, i):
            i_j = i_2 + j*j
            if i_j > N_2:
                break

            tuples += 1
            counts[i_j] += 8

    print("\ttuples:", tuples)
    print("\tsum(counts):", sum(counts))
    print("\tmax(counts):", max(counts))

    return counts


def get_n2_counts_factoring(N_2):
    import primesieve

    counts = array.array('H', [4]) * (N_2+1)
    assert counts.itemsize == 2, counts.itemsize
    assert len(counts) == N_2 + 1

    updates = 0
    counts[0] = 1

    # https://mathworld.wolfram.com/SumofSquaresFunction.html
    for p in primesieve.primes(2, N_2):
        if p == 2:
            # Special handling
            continue

        if p % 4 == 3:
            for exp in range(1, 100, 2):
                # zero any multiple not to an even power
                pp = p ** exp
                if pp > N_2:
                    break

                ppp = p * pp

                # Zero any odd power
                # any m * pp, m % p != 0
                for i_pp in range(0, N_2, ppp):
                    stop = min(i_pp + ppp, N_2 + 1)
                    for m in range(i_pp + pp, stop, pp):
                        counts[m] = 0
                        updates += 1


        if p % 4 == 1:
            # Need to handle each exponent
            for exp in range(1, 100):
                pp = p ** exp
                if pp > N_2:
                    break

                ppp = p * pp

                # any m * pp, m % p != 0
                for i_pp in range(0, N_2, ppp):
                    stop = min(i_pp + ppp, N_2 + 1)
                    for m in range(i_pp + pp, stop, pp):
                        counts[m] *= 1 + exp
                        updates += 1

    print("\tupdates:", updates)
    return counts


def get_n2_counts_on_demand_factoring(N_2):
    import primesieve

    T0 = time.time()

    factor = array.array('I', range(N_2+1))
    assert factor.itemsize == 4, factor.itemsize
    assert len(factor) == N_2 + 1

    r = math.isqrt(N_2)
    updates = 0
    for p in primesieve.primes(2, r+1):
          for m in range(p * p, N_2+1, p):
              factor[m] = p
              updates += 1

    # https://mathworld.wolfram.com/SumofSquaresFunction.html
    counts = array.array('H', [4]) * (N_2+1)
    assert counts.itemsize == 2, counts.itemsize
    assert len(counts) == N_2 + 1
    counts[0] = 1

    T1 = time.time()
    print("\tupdates:", updates)
    print(f"\tany factor done ({T1-T0:.2f} secs)")

    num_factors = 0
    for i, f in enumerate(tqdm(factor[2:]), 2):
        if f == 0:  # has factor of p
            continue

        factors = [f]
        n = i // f
        while n > 1:
            f = factor[n]
            n //= f
            num_factors += 1
            factors.append(f)

        # 4 * (1 + exp) for all exp for factors of form (4k+1)
        m = 4
        for p, e in Counter(factors).items():
            if p % 4 == 1:
                m *= 1 + e
            elif p % 4 == 3:
                if e % 2 == 1:
                    m = 0
        counts[i] = m

    T2 = time.time()
    print(f"\tsum of square representations done ({T2-T1:.2f} secs, total: {T2-T0:.2f})")
    print("\tnumber of factors:", num_factors)
    return counts


def enumerate_N2(N):
    counts1 = get_n2_counts(N)
    counts = counts1

    if N > 1000:
        # Nice verification check from A057655
        assert sum(counts[:1000 + 1]) == 3149

    if False:
        print()
        counts2 = get_n2_counts_on_demand_factoring(N)
        print()

        for i, (c1, c2) in enumerate(zip(counts1, counts2)):
            if c1 != c2:
                print ("mismatch at", i, "\t", c1, c2)

        assert(counts1 == counts2)


    # This code takes SO long but should be O(n) (AKA fast)
    PI = pi()

    A_n = 0
    record = 1 # To avoid initial zero term
    record_temp = record
    A000099 = []  # index
    A000036 = []  # diff
    A000323 = []  # A(index)
    for n, a in enumerate(counts):
        A_n += a

        V_n_temp = math.pi * n
        if abs(A_n - V_n_temp) + 0.01 < record_temp:
            # Avoid expensive decimal computation
            continue

        V_n = PI * n
        P_n = A_n - V_n
        # TODO fix this
        if abs(P_n) > record:
            nearest = P_n.to_integral()
            # Boring, as value is often pi
            improved = abs(P_n) - record
            A000099.append(n)
            A000036.append(nearest)
            A000323.append(A_n)
            nth = len(A000099)
            if (nth < 20) or (nth % 5 == 0) or (nth > 950):
                print(f"| {nth:3} | {n:11} | {nearest:14} | {A_n:14} |")
            record = abs(P_n)
            record_temp = float(record)

    for fn, An in [("b000099.txt", A000099), ("b000036.txt", A000036), ("b000323.txt", A000323)]:
        with open(fn, "w") as f:
            for i, a in enumerate(An, 1):
                f.write(f"{i} {a}\n")

# Memory usage is 2bytes per n, 2GB per 10^9
enumerate_N2(10736688777)  # 1000 records
