import math

from decimal import Decimal, getcontext
getcontext().prec = 50
getcontext().Emax = 999999999999999999

#import gmpy2

for n in range(1, 133+1):
    log2 = max(2**n - (n+1), 2**(n-1))
    log10 = Decimal(2).log10() * log2

    digits = int(log10 + 1)
    start = int(10 ** (6 + (log10%1)))
    END_MOD = 10 ** 6
    end1 = pow(2, 2**n - (n+1), END_MOD)
    end2 = (pow(2, n, END_MOD) - 1) * pow(2, 2**(n-1) - n, END_MOD)
    end = (end1 + end2) % END_MOD
    print (f"{n}   log2(an): {log2}   log10(an): {log10:.3f}   {start}...{end:06d}<{digits}>")

    if n <= 30:
        n = Decimal(n)
        #n = gmpy2.mpfr(n)
        an1 = 2**(2**n - (n+1))
        an2 = (2**n-1)*(2**(2**(n-1) - n))
        an = an1 + an2
        #print (an1 / an)
        print (n, an, "\t", an.log10())
