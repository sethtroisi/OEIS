import array
import bisect
import itertools
import math

from decimal import Decimal, getcontext
from tqdm import tqdm

# From https://docs.python.org/3/library/decimal.html
def pi():
    """Compute Pi to the current precision.

    >>> print(pi())
    3.141592653589793238462643383

    """
    getcontext().prec += 2  # extra digits for intermediate steps
    three = Decimal(3)      # substitute "three=3.0" for regular floats
    lasts, t, s, n, na, d, da = 0, three, 3, 1, 0, 0, 24
    while s != lasts:
        lasts = s
        n, na = n+na, na+8
        d, da = d+da, da+32
        t = (t * n) / d
        s += t
    getcontext().prec -= 2
    return +s               # unary plus applies the new precision


def get_n3_counts(N):
    counts = array.array('I', [0]) * (N+1)
    assert counts.itemsize == 4, counts.itemsize
    assert len(counts) == N + 1

    tuples = 0
    for i in tqdm(itertools.count(), total=math.isqrt(N) + 1):
        i_2 = i*i
        if i_2 > N:
            break

        for j in range(i+1):
            j_2 = j*j
            i_j = i_2 + j_2
            if i_j > N:
                break

            k_max = math.isqrt(N - i_j)
            assert k_max ** 2 + i_j <= N
            for k in range(min(k_max, j)+1):
                k_2 = k*k
                i_j_k = i_j + k_2
                #assert i_j_k <= N

                # This logic can be removed if j=i, k=i cases are handled
                # seperately, but it doesn't speed up code any
                same = (i == j == k) + (i == j) + (j == k)
                mult = (6, 3, None, 1)[same]

                zeros = (i == 0) + (j == 0) + (k == 0)
                sign = (8,4,2,1)[zeros]

                tuples += 1
                counts[i_j_k] += mult * sign

    print("\ttuples:", tuples)
    print("\tsum(counts):", sum(counts))
    print("\tmax(counts):", max(counts))

    return counts


def get_n3_counts_v2(N):
    """
    Get number of representations of n as (i*i + j*j + k*k) for all n <= N

    3-5x faster than v1
    """
    counts = array.array('I', [0]) * (N+1)
    assert counts.itemsize == 4, counts.itemsize
    assert len(counts) == N + 1

    r_2 = math.isqrt(N)

    tuples = 0

    # i == j == k == 0
    counts[0] += 1

    # i == j == k
    for i in range(1, math.isqrt(N // 3)+1):
        n = 3*i*i
        assert n <= N
        tuples += 1
        counts[n] += 8

    # i = j, j > k
    for i in range(1, math.isqrt(N // 2)+1):
        temp = 2*i*i
        assert temp <= N

        # k = 0
        tuples += 1
        counts[temp] += 3 * 4

        # k > 0, k < j
        for k in range(1, i):
            n = temp + k*k
            if n > N: break
            tuples += 1
            counts[n] += 24  # 3 * 8

    # i > j = k
    for i in range(1, r_2+1):
        i_2 = i*i
        assert i_2 <= N

        # j = k = 0
        tuples += 1
        counts[i_2] += 6  # 3 * 2

        # j = k, j > 0
        for j in range(1, i):
            n = i_2 + 2*j*j
            if n > N: break
            tuples += 1
            counts[n] += 24  # 3 * 8

    for i in range(1, r_2 + 1):
        i_2 = i*i
        # i > j, k = 0
        for j in range(1, i):
            j_2 = j*j
            i_j = i_2 + j_2
            if i_j > N:
                break
            # k = 0
            tuples += 1
            counts[i_j] += 24  # 6 * 4

    # Roughly 9/10 N upper bound 5/4 N
    print(f"\tprocessed {tuples:,} special case pairs before main loop")

    # Can build sorted list of (j, k) pairs to better help with locality
    pairs = [] #array.array('I')
    for i in tqdm(range(1, r_2 + 1)):
        i_2 = i*i

        # Remove anything where i_2 + pair > N
        if pairs and pairs[-1] + i_2 > N:
          del pairs[bisect.bisect(pairs, N - i_2):]
          assert not pairs or pairs[-1] + i_2 <= N

        # add new pairs for j = i-1
        j = i - 1
        j_2 = j*j
        i_j = i_2 + j_2
        if i_j + 1 <= N:
          merge_start = bisect.bisect_left(pairs, j_2)

          for k in range(1, j):
              k_2 = k*k
              n = i_j + k_2
              if n > N: break
              pairs.append(j_2 + k_2)

          #pairs[merge_start:] = sorted(pairs[merge_start:])
          pairs.sort()

          # Merge should be O(n) but is slower in practice
          #pairs = list(heapq.merge(pairs, new_pairs))
          #assert pairs == sorted(pairs), (pairs[:merge_start], pairs[merge_start:], pairs)

        # TODO: check if duplicate items ever appear in pairs
        # TODO: try using bitarray & bitscan.

        # A clever way to get progress reports from tqdm
        if i > 0 and ((i * 20) % r_2) < 20:
          print(f"\t{i=:<6d}  |pairs|={len(pairs)}, processed={tuples:,}")

        tuples += len(pairs)

        # Data being sorted means access to counts is sequential and fast
        for p in pairs:
            n = i_2 + p
            assert n <= N
            counts[n] += 48  # 6 * 8

    print("\ttuples:", tuples)
    print("\tsum(counts):", sum(counts))
    print("\tmax(counts):", max(counts))

    return counts


def enumerate_n3(N):
    counts = get_n3_counts_v2(N)
    #assert counts == get_n3_counts(N)

    if N > 1000:
        # Nice verification check from A117609
        sum_1000 = sum(counts[:1000 + 1])
        assert sum_1000 == 132451, sum_1000

    getcontext().prec = int(2 * math.log10(N) + 10)
    M = Decimal(4) / 3 * pi()

    def V(n):
        """Slightly higher precision"""
        return M * Decimal(n ** 3).sqrt()

    print(f'| {"nth":3} | {"n = A000092":11} | {"P(n) = A000223":14} | {"A(n) = A000413":14} |')

    A_n = 0
    record = 1 # To avoid initial zero term
    record_float = float(record)
    A000092 = []
    A000223 = []
    A000413 = []
    for n, a in enumerate(counts):
        A_n += a

        V_n_float = 4/3 * math.pi * n ** (3/2)
        P_n_temp = abs(A_n - V_n_float)
        if P_n_temp + 1 < record_float:
            # Max error is floating point imprecision of V_n.
            # No need to invoke expensive V(n)
            continue

        V_n = V(n)
        P_n = (A_n - V_n)
        if abs(P_n) > record:
            # A000223 uses rounded version of V(n)
            # Check if floating point rounds differently
            if V_n.to_integral() != round(V_n_float):
                print(f"float based V_n would round differently at {n}  {V_n} vs {V_n_float}")

            # Calculate rounded version of P(n)
            P_n_rounded = (A_n - V_n).to_integral()

            A000092.append(n)
            A000223.append(P_n_rounded)
            A000413.append(A_n)
            record = abs(P_n)
            record_float = float(P_n)
            nth = len(A000092)
            if (nth < 10) or (nth % 5 == 0) or (nth in range(50,55)) or (nth > 120):
                print(f"| {nth:3} | {n:11} | {P_n_rounded:14} | {A_n:14} |")

    # If more records than on file
    if len(A000092) > 188:
        for fn, An in [("b000092.txt", A000092), ("b000223.txt", A000223), ("b000413.txt", A000413)]:
            with open(fn, "w") as f:
                for i, a in enumerate(An, 1):
                    f.write(f"{i} {a}\n")



# 100 terms in 1 second
#enumerate_n3(1560000)
# 124 terms in 11 seconds
enumerate_n3(10 * 10 ** 6)

# 131 terms in 27 seconds
#enumerate_n3(20 * 10 ** 6)

# 151 terms in 166 seconds
#enumerate_n3(63 * 10 ** 6)

# 188 terms in <45 minutes
#enumerate_n3(450 * 10 ** 6)

# 194 in <76 minutes with pypy3
#enumerate_n3(600 * 10 ** 6)

# 210 in <275 minutes with pypy3
#enumerate_n3(1400 * 10 ** 6)
