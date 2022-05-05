# OEIS Quadratic Forms Tracker

I (Seth) have worked on a few of quadratic representation/population sequences namely.

* [A000047](https://oeis.org/A000047) - `x^2 - 2*y^2` - prime signature, p in (3,5) mod 8 occur to even power
* [A000049](https://oeis.org/A000049) - `3*x^2 + 4*y^2` - ???
* [A000050](https://oeis.org/A000050) - `x^2 + y^2` - prime signature, p = 3 mod 4 occur to even power
* [A000067](https://oeis.org/A000067) - `x^2 + 2*y^2` - prime signature,  p in (5,7) mod 8 occur to even power
* [A000072](https://oeis.org/A000072) - `x^2 + 4*y^2` -
    * A000072 can be produced from A000050 = `B_4` via `B_4(n) = B_1(n) - B_1(n-1) + B_1(n-2)`
* [A000074](https://oeis.org/A000074) - odds `x^2 + y^2` - first derivative of A000050!
* [A000205](https://oeis.org/A000205) - `x^2 + 3*y^2` - prime signature, p = 2 mod 3 occur to even power

This covers `B_n()` with n = `{-2, 1, 2, 3, 4}`.

* Population vs Enumeration
    * [A000047](https://oeis.org/A000047) is population of [A035251](https://oeis.org/A035251)
    * [A000049](https://oeis.org/A000049) is population of [A020677](https://oeis.org/A020677)
    * [A000050](https://oeis.org/A000050) is population of [A001481](https://oeis.org/A001481)
    * [A000067](https://oeis.org/A000067) is population of [A002479](https://oeis.org/A002479)
    * [A000205](https://oeis.org/A000205) is population of [A003136](https://oeis.org/A003136)

* No Population but easy to count
    * [A031363](https://oeis.org/A031363) - `5*x^2 - y^2` - has a prime signature p in (2,3) mod 5

## TODO

* Class No. 2 Sequences
    * Determine what primes appear in P0 / P1
        * For `n=6` this is `P0 = [p for p in primes if (2/p) == 1 and (-3/p) == 1]`
            * Will need to figure out how to init counts (and make sure they work the same)
    * Need to figure out how to handle multiplicity of P1
        * Maybe same exact way as inclusion - exclusion?

* Add sequence for population of [A020670](https://oeis.org/A020670) - `x^2 + 7*y^2` (has prime representation)
    * Let R(n) = [2][P + 7][Q^2]
    * A(n) = R(n) - R(N-1) + R(N-2) to handle 2 not being represented


---

Sequences to consider

* [A000018](https://oeis.org/A000018) - `x^2 + 16*y^2`
* [A000021](https://oeis.org/A000021) - `x^2 + 12*y^2`
* [A000024](https://oeis.org/A000024) - `x^2 + 10*y^2` - `B_10(n)` is promised in Class No. 2
* [A000077](https://oeis.org/A000077) - `x^2 + 6*y^2`  - `B_6(n)` and `B_2,3(n)` might be solvable at the same time? see `(38)`
* [A054150](https://oeis.org/A054150) - `x^2 + 5*y^2` - `B_5(n)` is promised in Class No. 2
* [A000075](https://oeis.org/A000075) - `2*x^2 + 3*y^2` - `B_2,3(n)` See A000077 above
* [A000286](https://oeis.org/A000286) - `2*x^2 + 5*y^2` - `B_2,5(n)` is promised in Class No. 2
* [A000049](https://oeis.org/A000049) - `3*x^2 + 4*y^2`
    * `B_3,4(n)` can be produced from `B_3(n)` (known) and `B_12(n)` (A21)
    * `E_3,4(n)` = `B_3(n-2)` (known)
* [A000076](https://oeis.org/A000076) - `4*x^2 + 4*x*y + 5*y^2`
    * `B_4,4,5(n)` can be derived at the same time as `B_16(n)` and `B1(n)`

---

Reading D. Shanks and L. P. Schmid, [Variations on a theorem of Landau](
http://dx.doi.org/10.1090/S0025-5718-1966-0210678-1). Part I, Math. Comp., 20
(1966), 551-569. it seems likely all of the small number sequences have clever
prime representations (or recurrances).

Sequences 72, 75, 76, and 77 are all populations (# of terms < 2^n) of
[quadratic forms](https://oeis.org/wiki/Index_to_OEIS:_Section_Qua#quadpop).
None have known (to me) prime representations. So the best way I know is
`O(n^2)` enumeration of all matching points.




## Solved sequences

### A000047 [Project README](../A000047/README.md)

Quoting [A035251](https://oeis.org/A035251)

> A positive number n is representable in the form x^2 - 2y^2 iff every prime p == 3 or 5 (mod 8) dividing n occurs to an even power.

### A000050 [Project README](../A000050/README.md)

From P. Moree and H. J. J. te Riele,
[The hexagonal versus the square lattice, arXiv:math/0204332](https://arxiv.org/abs/math/0204332)
[math.NT], 2002. I learned

> Lemma 1. A positive integer n is represented by the form X^2 + Y^2 if and only if
every prime factor p of n of the form p = 3 (mod 4) occurs to an even power.
> A positive integer n is represented by the form X^2 + 3\*Y^2 if and only if every prime
factor p of n of the form p = 2 (mod 3) occurs to an even power.


Extras

[Kaplansky's theorem on quadratic forms (wikipedia)]
(https://en.wikipedia.org/wiki/Kaplansky%27s_theorem_on_quadratic_forms)

David Brink, [Five peculiar theorems on simultaneous representation of primes by quadratic
forms, doi:10.1016/j.jnt.2008.04.007](https://doi.org/10.1016/j.jnt.2008.04.007)

* 4k + 1 even multiplicity (A268377) is almost the same as A054165 (3 x^2 + 8 y^2)
  * See A020680

---

