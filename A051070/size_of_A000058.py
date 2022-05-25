from decimal import Decimal, getcontext


getcontext().Emax=2**58

an = Decimal(2)
an_int = 2
for n in range(1, 58+1):
    an = an ** 2 - an + 1
    an_int = an_int ** 2 - an_int + 1
    if n >= 5:
        an_int %= 10 ** 6
    print (n, an, "...", an_int)
