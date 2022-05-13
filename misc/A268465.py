import gmpy2
import primesieve

it = primesieve.Iterator()

a, b, c = 0, 0, 0
An = []
m = 1
while True:
    a, b, c = b, c, m * it.next_prime()
    m += 1
    if m <= 3:
        continue

    n = a + b + c
    if gmpy2.is_prime(n):
        An.append(n)
        if len(An) == 100:
            break

with open("b268465.txt", "w") as f:
    for n, an in enumerate(An, 1):
        f.write(f"{n} {an}\n")
        print(n, an)
