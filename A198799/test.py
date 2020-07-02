import sympy

mod1 = [p for p in sympy.primerange(2,50000) if p % 6 == 1]
print (mod1[:10])

# A198799 (smallest m that can be written in EXACTLY n ways)
# These seem to return 2*m, except sometimes 2*m-1 where one solution seems to be "lost" because x=0?
A198799 = "49 637 1729 8281 12103 1529437 53599 157339 593047 19882681 375193 68574961 2989441 7709611 1983163 47738317081 4877509 21169376772835837 18384457 377770939"

# A198774 (3 reps), All of these return m=6
A198774 = """1 637
2 931
3 1183
4 1519
10 2548
20 4693
50 9751
100 18228
200 35017
400 69727
600 105889
800 142492
1000 179977"""

# A198775 (4 reps), All of these return m=8
A198775 = """1 1729
2 2821
3 3367
94 22204
106 24583
115 25753
128 27937
142 30121
158 32116
173 34333
193 37387
211 39711
224 41587
239 43491
249 44499
250 44548"""

A = A198799
A = list(map(int, A.split(" ")))

#A = A198775
#A = [int(line.split(" ")[1]) for line in A.split("\n")]

print ()
for i, a in enumerate(A, 2):
    factors = sympy.factorint(a)
    #assert all((p in mod1) or (p == 3) for p in factors), factors

    exponents = list(c for p, c in factors.items() if p % 6 == 1)
    m = sympy.prod([e+1 for e in exponents])

    print (f"{i}\t{a:18} {str(exponents):20} {m:2} {str(factors)}")


def brute_verify(k):
    M = int(2 * (k/3) ** 0.5 + 5)
    count = 0
    for y in range(M):
        for x in range(y):
            if x*x + y*y + x*y == k:
                count += 1
                print (count, x, y)
    print()
    return count

print()
brute_verify(49)
brute_verify(637)
brute_verify(1729)
brute_verify(8281)
brute_verify(12103)
