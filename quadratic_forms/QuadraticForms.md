# OEIS Quadratic Forms Tracker

I (Seth) have worked on a few of quadratic representation/population sequences namely.

* [A000047](https://oeis.org/A000047) - `x^2 - 2*y^2`
* [A000049](https://oeis.org/A000049) - `3*x^2 + 4*y^2`
* [A000050](https://oeis.org/A000050) - `x^2 + y^2`

## A000047 [Project README](../A000047/README.md)

Quoting [A035251](https://oeis.org/A035251)

> A positive number n is representable in the form x^2 - 2y^2 iff every prime p == 3 or 5 (mod 8) dividing n occurs to an even power.

## A000049 [Project README](../A000049/README.md)

OEIS sequence for representation [A020677](https://oeis.org/A020677)

Working on signature in [quadratic\_representations.py](quadratic_representations.py)

* `3*x^2 + 4*y^2` with even `x`, `e.g.` `x = 2*z` -> `3*4*z^2 + 4*y^2 = 4*(3*z^2 + y^2)`
  * See A000050 -> 4 * A000050

* with `p` to even power can come from `(x/p, y/p)`

* when `z = 3*x^2 + 4*y^2 = p*q` things make sense
  * wlog (p == 3 or p = 7 mod 12) and (q = 1 mod 12)

* Things break down when `z = p^3 * q`

## A000050 [Project README](../A000050/README.md)

From P. Moree and H. J. J. te Riele,
[The hexagonal versus the square lattice, arXiv:math/0204332](https://arxiv.org/abs/math/0204332)
[math.NT], 2002. I learned

> Lemma 1. A positive integer n is represented by the form X^2 + Y^2 if and only if
every prime factor p of n of the form p = 3 (mod 4) occurs to an even power.
> A positive integer n is represented by the form X^2 + 3\*Y^2 if and only if every prime
factor p of n of the form p = 2 (mod 3) occurs to an even power.
