import sympy
import primesieve

it = primesieve.Iterator()
n = 0
while True:
    np = it.next_prime() + 1
    factors = sympy.factorint(np)
    if len(factors) == 5:
        n += 1
        print(n, np)
        #if np > 17000:
        #    break
