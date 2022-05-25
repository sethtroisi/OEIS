import math
from decimal import Decimal, getcontext
import itertools
import gmpy2

def cube_root(A):
    """Find the cube root via newton's method"""
    x0 = (A-1)/3
    xn = (2 * x0 + A / (x0*x0) ) / 3

    # how to check if error is in last decimal place?
    for iteration in range(1000):
        if xn == x0:
            break
        x0 = xn
        xn = (2 * x0 + A / (x0*x0) ) / 3
    #print (A, xn - x0)
    return xn


record = 0.01

def A000148(n):
    global record
    max_xi = math.isqrt(n ** 3)
    #print (f"\tn: {n}, max(x_i): {max_xi}")
    assert max_xi ** (2/3) <= n

    # Split into ingeters and irrationals
    rational = []
    irrational = list(range(1, max_xi+1))
    for i in range(1, max_xi+1):
        if i ** 2 > n: break
        rational.append(i ** 2)
        irrational.remove(i ** 3)

    irrational = [CUBE_ROOTS_OF_SQUARES[i] for i in irrational]

    count = 0
    for sum_two in map(sum, itertools.combinations_with_replacement(rational, r=2)):
        if sum_two <= n:
            count += 1

    ii = len(irrational)-1  # irrational index
    for r in rational:
        while ii >= 0 and irrational[ii] > (n - r):
            ii -= 1
        if ii < 0:
            break
        count += (ii + 1)

    # TODO technically could use continued fractions
    bi = len(irrational)-1 # irrational index
    for ai, a in enumerate(irrational):
        test = n - a
        if test <= a:
            break

        error = abs(irrational[bi] - test)
        if error < record:
            record = error
            t1 = CUBE_ROOTS_OF_SQUARES.index(a)
            t2 = CUBE_ROOTS_OF_SQUARES.index(irrational[bi])
            print(f"\t{error:.2g} from {t1}^(2/3) + {t2}^(2/3)")

        while bi >= ai and irrational[bi] > test:
            bi -= 1

        if bi < ai:
            break

        count += bi - ai + 1

    return count



MAX_XI = math.isqrt(1000 ** 3)

getcontext().prec = 60
CUBE_ROOTS_OF_SQUARES = [Decimal(0), Decimal(1)] + [cube_root(Decimal(i ** 2)) for i in range(2, MAX_XI+1)]
#CUBE_ROOTS_OF_SQUARES = [gmpy2.mpfr(i) ** (2/3) for i in range(MAX_XI+1)]


with open("b000148.txt", "w") as bfile:
    for n in range(2, 1000+1):
        an = A000148(n)
        bfile.write(f"{n} {an}\n")
        print(n, an)
