#!/usr/bin/env python

import sys


def count_quadratic_form(n:int, a: int, b: int, c: int) -> int:
    """count of all numbers n = a x^2 + b x*y + c y^2 <= N."""
    population = set()

    assert c >= 0, (a, b, c)

    for x in range(n):
        temp_x = a * x*x
        if temp_x > n:
            break

        for y in range(n):
            temp = temp_x + b * x*y + c * y*y
            if temp > n:
                break

            if temp > 0:
                population.add(temp)

    population.discard(0)

    if n == 2 ** 24:
        singles = set()
        import sympy
        for p in sorted(population):
            fs = sympy.factorint(p)
            print("\t", p, "\t", fs)
            for f, c in fs.items():
                if c in (2,4):
                    singles.add(f)

        print(sorted(singles))


    return len(population)


def main():
    args = list(map(int, sys.argv[1:]))

    a, b, c = 0, 0, 0

    if len(args) == 1:
        a, c = 1, args[0]
        name = f"x^2 + {c} y^2"
    elif len(args) == 2:
        a, c = args
        name = f"{a} x^2 + {c} y^2"
    elif len(args) == 3:
        a, b, c = args
        name = f"{a} x^2 + {b} x*y + {c} y^2"

    print(f"Enumerating  {name} <= 2^n")
    for n in range(1, 20):
        count = count_quadratic_form(2 ** n, a, b, c)
        print(f"| {n:2d} | {name} | {count:8} |")



if __name__ == "__main__":
    main()
