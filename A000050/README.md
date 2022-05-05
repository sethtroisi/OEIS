# [A000050](https://oeis.org/A000050)

From P. Moree and H. J. J. te Riele,
[The hexagonal versus the square lattice, arXiv:math/0204332](https://arxiv.org/abs/math/0204332)
[math.NT], 2002. I learned

> Lemma 1. A positive integer n is represented by the form X^2 + Y^2 if and only if
every prime factor p of n of the form p = 3 (mod 4) occurs to an even power. A
positive integer n is represented by the form X^2 + 3\*Y^2 if and only if every prime
factor p of n of the form p = 2 (mod 3) occurs to an even power.

This allows for a simple adaption of the prime signature based code.

[A000072](https://oeis.org/A000072) and [A000074](https://oeis.org/A000074) are
both derived from this sequence with

```
# For writing out b000050.txt
$ cat README.md  | awk '/\| .. \| [1-9]/ {print $2, $4}'

# A72(n) = A50(n) - A50(n-1) + A50(n-2)
$ cat b000050.txt | awk '{if ($1) print($1, $2-l+ll); ll=l; l=$2;}' | tee b000072.txt

# A74(n) = A50(n) - A50(n-1)
$ cat b000050.txt | awk '{if ($1) print($1, $2-l); l =$2;}' | tee b000074.txt
```

This sequence only had 35 terms extending it to 50 unblocks various OEIS
self referential sequences.

See [A051070](https://oeis.org/A051070) and
[Self-referential sequences (wikipedia)](https://en.wikipedia.org/wiki/On-Line_Encyclopedia_of_Integer_Sequences#Self-referential_sequences)

Initially I adapted (nearly verbatim) `A000050_queue` and `A000050_segmented_hash`
from [OEIS/A000049](../A000049).


## Results

| n  | a(n)            | iters          | queue time(s)| hash (s) |
|----|-----------------|----------------|--------------|----------|
| 25 | 6402706         | 13181716       | 0.8     secs |          |
| 26 | 12534812        | 26360545       | 1.7     secs |          |
| 27 | 24561934        | 52717056       | 3.6     secs |          |
| 28 | 48168461        | 105428281      | 7.5     secs |          |
| 29 | 94534626        | 210848467      | 15.5    secs |          |
| 30 | 185661958       | 421685278      | 31.8    secs | 0.9      |
| 31 | 364869032       | 843354269      | 65.9    secs | 1.5      |
| 32 | 717484560       | 1686685460     | 138.0   secs | 3.3      |
| 33 | 1411667114      | 3373338369     | 290.5   secs | 7.1      |
| 34 | 2778945873      | 6746630485     | 618.5   secs | 14.4     |
| 35 | 5473203125      | 13493195702    |              | 27.2     |
| 36 | 10784591780     | 26986298904    |              | 49.4     |

| n  | a(n)               | count primes = 3 mod 4 <= 2^n | time (s) |
|----|--------------------|------------------|--------------|----------|
| 37 | 21259581728        | 2793255158       | 0.78     |
| 38 | 41926021639        | 5433148653       | 1.30     |
| 39 | 82714309017        | 10575966554      | 2.18     |
| 40 | 163243760848       | 20601565111      | 3.59     |
| 41 | 322287239923       | 40158327450      | 5.84     |
| 42 | 636491338440       | 78330571530      | 9.86     |
| 43 | 1257412026877      | 152880882090     | 16.14    |
| 44 | 2484802305138      | 298558278909     | 26.84    |
| 45 | 4911668621134      | 583373452992     | 44.94    |
| 46 | 9711438405328      | 1140499484946    | 76.54    |
| 47 | 19206579700090     | 2230816749376    | 127.55   |
| 48 | 37994743870997     | 4365594657172    | 206.51   |
| 49 | 75179493114567     | 8547216751556    | 346.41   |
| 50 | 148789852637074    | 16741690493299   | 584.42   |
| 51 | 294537921458811    | 32806450071861   | 972.88   |
| 52 | 583175525599178    | 64312752341407   | 1626.29  |
| 53 | 1154898832332412   | 126126353622230  | 2737.53  |
| 54 | 2287556346291550   | 247445104415883  | 4589.74  |
| 55 | 4531893363844648   | 485634974839804  | 7770.19  |
| 56 | 8979757047285127   | 953439696436532  | 13321.32 |
| 57 | 17796052510039673  | 1872505595140426 | 22467.86 |
| 58 | 35273964779995151  | 3678700139357953 (18652 secs) | 37261.30 |
| 59 | 69928440727755836  | 7229396458702824 (31166 secs) | 63353.79 |
| 60 | 138650032775784383 | 14211547255716680 | |

### Dreaming of A(72) and A(74)

`a(60)` used 53 GB of memory vs the theoretical 48 GB `3 * int64 * 2 * (2^60)^0.5 / (8 * 1024^3)`.
With a brief peak to `4/3 * 53 =` 70 GB.

Memory usage can be decreased by 1/3 by removing `v` values from `counts\_backing`.
And peak usage can be reduced by using `pop_back()` and `shrink_to_fit`

If both were implemented `A72(72)` could be computed via `A50(72)` with `48 * 2^(72/2 - 30)` =
3072 GB of RAM! An amount easily accessible in cloud based machines.

Runtime is `O((2^n)^(3/4))`, extrapolating from `a(55)` ` yields an estimate of 109 CPU-years.
Replacing the PrimePi code with [primecount](https://github.com/kimwalisch/primecount) could reduce
runtime by 30-50% but requires some work to figure out how to handle the two classes
`p mod 3 = 2 vs 3` and to extract all intermediate values (primepi(n/2), primepi(n/3).... In a
private communication Kim warned this would be a sizable undertaking.


```
# Fast
$ g++ -O3 -std=c++17 A000050_signature.cpp ../utils/count_special_primes.cpp -lprimesieve

# For double check / historical interest
$ g++ -O3 A000050_queue.cpp
$ g++ -O3 -march=native -fopenmp -Wall -Werror -std=c++17 A000050_segmented_hash.cpp
```

