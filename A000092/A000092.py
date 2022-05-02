import array
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


def get_n3_counts(N_2):
    # I wish array had something like numpy.zeros

    counts = array.array('I', [0] * (N_2+1))
    assert counts.itemsize == 4, counts.itemsize
    tuples = 0

    # 1.5x slower
    #counts = array.array('I')
    #counts.frombytes(bytearray(4 * (N_2+1)))

    # Why is this 2-3x slower?
    #counts = array.array('I', (0 for i in range(N_2+1)))

    assert len(counts) == N_2 + 1

    for i in tqdm(itertools.count(), total=math.isqrt(N_2) + 1):
        i_2 = i*i
        if i_2 > N_2:
            break

        for j in range(i+1):
            j_2 = j*j
            i_j = i_2 + j_2
            if i_j > N_2:
                break

            for k in range(j+1):
                k_2 = k*k
                i_j_k = i_j + k_2
                if i_j_k > N_2:
                    break

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


def enumerate_n3(N_2):
    counts = get_n3_counts(N_2)

    if N_2 > 1000:
        # Nice verification check from A117609
        assert sum(counts[:1000 + 1]) == 132451

    getcontext().prec = int(2 * math.log10(N_2) + 10)
    M = Decimal(4) / 3 * pi()

    def V(n):
        """Slightly higher precision"""
        return M * Decimal(n ** 3).sqrt()

    print(f'| {"nth":3} | {"n = A000092":11} | {"P(n) = A000223":14} | {"A(n) = A000413":14} |')

    A_n = 0
    record = 1 # To avoid initial zero term
    A000092 = []
    A000223 = []
    A000413 = []
    for n, a in enumerate(counts):
        A_n += a

        V_n_float = round(4/3 * math.pi * n ** (3/2))
        P_n_temp = abs(A_n - V_n_float)
        if P_n_temp + 10 < record:
            # No need to invoke expensive V(n)
            continue

        V_n = V(n).to_integral()
        if V_n != V_n_float:
            print(f"Mismatch V_n at {n}  {V_n} vs {V_n_float}")

        P_n = abs(A_n - V_n)
        if P_n > record:
            A000092.append(n)
            A000223.append(P_n)
            A000413.append(A_n)
            record = P_n
            nth = len(A000092)
            if (nth < 20) or (nth % 5 == 0) or (nth > 90):
                print(f"| {nth:3} | {n:11} | {P_n:14} | {A_n:14} |")

    for fn, An in [("b000092.txt", A000092), ("b000223.txt", A000223), ("b000413.txt", A000413)]:
        with open(fn, "w") as f:
            for i, a in enumerate(An, 1):
                f.write(f"{i} {a}\n")



# For 100 terms in 1 second
#enumerate_n3(1560000)

enumerate_n3(2 * 10 ** 7)
