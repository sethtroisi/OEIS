# OEIS Quadratic Forms Tracker

I (Seth) have worked on a few of quadratic representation/population sequences namely.

* [A000047](https://oeis.org/A000047) - `x^2 - 2*y^2` - prime signature, p in (3,5) mod 8 occur to even power
* [A000049](https://oeis.org/A000049) - `3*x^2 + 4*y^2` - ???
* [A000050](https://oeis.org/A000050) - `x^2 + y^2` - prime signature, p = 3 mod 4 occur to even power
* [A000074](https://oeis.org/A000074) - odds `x^2 + y^2` - first derivative of A000050!
* [A000205](https://oeis.org/A000205) - `x^2 + 3*y^2` - prime signature, p = 2 mod 3 occur to even power

TODO extend 205 a few terms

---

Sequences to consider

* [A000018](https://oeis.org/A000018) - `x^2 + 16*y^2`
* [A000021](https://oeis.org/A000021) - `x^2 + 12*y^2`
* [A000067](https://oeis.org/A000067) - `x^2 + 2*y^2` - has a prime signature p in (5,7) mod 8
* [A000072](https://oeis.org/A000072) - `x^2 + 4*y^2`
* [A000075](https://oeis.org/A000075) - `2*x^2 + 3*y^2`
* [A000076](https://oeis.org/A000076) - `4*x^2 + 4*x*y + 5*y^2`
* [A000077](https://oeis.org/A000077) - `x^2 + 6*y^2`
* [A000286](https://oeis.org/A000286) - `2*x^2 + 5*y^2`
* [A031363](https://oeis.org/A031363) - `5*x^2 - y^2` - has a prime signature p in (2,3) mod 5

## A000047 [Project README](../A000047/README.md)

Quoting [A035251](https://oeis.org/A035251)

> A positive number n is representable in the form x^2 - 2y^2 iff every prime p == 3 or 5 (mod 8) dividing n occurs to an even power.

## A000049 [Project README](../A000049/README.md)

OEIS sequence for representation [A020677](https://oeis.org/A020677)

Working on signature in [quadratic\_representations.py](quadratic_representations.py)

* `3*x^2 + 4*y^2` with even `x`, e.g. `x = 2*z` -> `3*4*z^2 + 4*y^2 = 4*(3*z^2 + y^2)`
  * 4 * any term in [A001481](https://oeis.org/A001481)

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
