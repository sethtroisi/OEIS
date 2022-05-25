# [A051070](https://oeis.org/A051070) "a(n) is the n-th term in sequence A_n"

See also [Self-referential sequences (Wikipedia)](
https://en.wikipedia.org/wiki/On-Line_Encyclopedia_of_Integer_Sequences#Self-referential_sequences)


### Background

I believe I found this sequence via a Project Euler problem related to
[A003136](https://oeis.org/A003136) or [A088534](https://oeis.org/A088534).

I had previously seen the OEIS self referencing sequence and at some point I looked it up and found
it was blocked at the 47th term (e.g. [A000047](https://oeis.org/A000047). From there I repurposed
some code I had written in college to enumerate taxicab numbers ([OEIS/A001235](../A001235)) to
enumerate initial terms of A000047. This spiralled out of control taking up two weeks of my life
while I obsessively optimized my enumeration of various quadratic forms.


### Updates

I have unblocked these sequences

1. [A000047](https://oeis.org/A000047) - Added terms 38-53
1. [A000049](https://oeis.org/A000049) - Added terms 37-49
    * This was a tought fight, My initial estimate was **150 days** for my taxicab
      (`priority_queue` based) code to finish. This was actually an **underestimate** with
      how enumeration slows down at larger N.
    * I found two clever tricks to dramatically speed this up
        1. Using `bitset<...>` as a low constant `O(1)` hash set.
        2. Breaking enumeration up into disjoint congruence classes.
1. [A000050](https://oeis.org/A000050) - Added terms 36-50
1. [A000092](https://oeis.org/A000092) - Added terms 49-131

Partially unblocked

1. [A000067](https://oeis.org/A000067) - Partially unblocked by calculating terms 36-55
    * Current code would take ~70 days to find `A000067(67)` and 100+ days for all intermediate terms
    * there's an possible 3x speedup from improving `prime_pi`
        * or using [primecount](https://github.com/kimwalisch/primecount) directly
    * Code can probably be multi-threaded (4-6x speedup)
1. [A000072](https://oeis.org/A000072) - Same as [A000067](https://oeis.org/A000067) see above.
1. [A000074](https://oeis.org/A000074) - Same as [A000067](https://oeis.org/A000067) see above.

https://oeis.org/A000072


### Blocking terms

1. [A000053](https://oeis.org/A000053) - finite (with less than 53 terms) -> -1
1. [A000054](https://oeis.org/A000054) - finite (with less than 53 terms) -> -1
1. [A000058](https://oeis.org/A000058) - `a(13)` has 1668 digits, `a(58)` will have ~ 58669977298272600 digits!
    * I think there's an argument for this number being "knowable" but write-down-able
1. [A000066](https://oeis.org/A000066) - Only 12 terms are know and problem class is "hard" (graph enumeration)
    * Last progress was 2012. Probably roughly `O(n!)`.
    Unlikely we'll find `a(20)` let alone `a(66)`
1. [A000067](https://oeis.org/A000067) - 50 known terms. [^1] - XX cpu-days
1. [A000072](https://oeis.org/A000072) - 35 known terms. [^1] - Theoretically possible (e.g. X cpu-years)
1. [A000074](https://oeis.org/A000074) - 50 known terms. [^1] - Theoretically possible (e.g. XX cpu-years)
1. [A000075](https://oeis.org/A000075) - 35 known terms. [^2] - XXX cpu-years, lots of ram, X weeks of coding
1. [A000076](https://oeis.org/A000076) - 39 known terms. [^2] - XXX cpu-years
1. [A000077](https://oeis.org/A000077) - 35 known terms. [^2] - XXX cpu-years
1. [A000080](https://oeis.org/A000080) - 25 known terms. Current code takes roughly `O(3.4^n)`
    * Does include working code (from 2015) with descriptive notes, thanks Martin Fuller!
1. [A000088](https://oeis.org/A000088) - 88 terms listed (but 0 to 87). `a(87)` is know (but has 1019 digits).
1. [A000101](https://oeis.org/A000101) - 80 terms.
    * I spent a lot of time working on this as part of a [GIMPS subproject](
https://www.mersenneforum.org/forumdisplay.php?f=131)

Large but easily computable terms

1. A00058(58) - 192523...920807<58669977298272603>
1. A000133(133) - 1200227...406272<3277926382039659640000387496003303446363>


[^1]: Has existing code that runs in `O((2^n)^(3/4))`, could be unblocked with more CPU time.

[^2]: Naive current code is `O(2^n)`. Implementing Class No. 2 and Section 9 from "Variations
      on a theorem of Landua" would likely improve this to `O((2^n)^(3/4))`. See
      [QuadraticForms.md](../quadratic_forms/QuadraticForms.md).
