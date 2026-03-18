# [A000047](https://oeis.org/A000047)

Quoting [A035251](https://oeis.org/A035251)

> A positive number n is representable in the form x^2 - 2y^2 iff every prime p == 3 or 5 (mod 8) dividing n occurs to an even power.

When I started this sequence only had 37 terms, extending it to 50 unblocks various OEIS self referential sequences.

See [A051070](https://oeis.org/A051070) and
[Self-referential sequences (wikipedia)](https://en.wikipedia.org/wiki/On-Line_Encyclopedia_of_Integer_Sequences#Self-referential_sequences)

I've adapted Lucy\_Hedgehog's [Count Primes](https://math.stackexchange.com/a/2283829/87805) to count the number of primes p == 3 or 5 (mod 8) under `n, n/2, n/3, ... ,n/sqrt(n), ... 5, 4, 3, 2, 1`  then I use inclusion-exclusion to count how many have the correct type of signature.

* 80-95% of the time is spent in `get_special_prime_counts` which I don't see how to parellize. So this will remain single threaded.

## Results

|n  | A000047(n)       | count primes     |time (c++)| python |
|---|------------------|------------------|----|---|
|38 | 37694785327      | 5433137518       | 1s | 23s (python) |
|39 | 74373523152      | 10575969498      | 1s | 34s (python) |
|40 | 146795295349     | 20601560352      | 2s | 62s (python) |
|41 | 289837406442     | 40158311317      | 2s | 87s (python) |
|42 | 572450196691     | 78330556749      | 4s | |
|43 | 1130980093866    | 152880952819     | 6s | |
|44 | 2235114244980    | 298558223205     | 14s | |
|45 | 4418409686942    | 583373483173     | 17s | |
|46 | 8736714417829    | 1140499468654    | 38s | 23m (python) |
|47 | 17279890334179   | 2230816679254    | 51s | 33m (python) |
|48 | 34185319295647   | 4365594817954    | 89s | 65m (python) |
|49 | 67645602070024   | 8547216668701    | 135s | 132m (python) |
|50 | 133886409875486  | 16741689949917   | 4m  194m (python) |
|51 | 265049159666520  | 32806450285778   | 6m  | |
|52 | 524814250058718  | 64312752938850   | 12m | |
|53 | 1039370837321432 | 126126352240561  | 17m | |
|54 | 2058817412464188 | 247445105913944  | 35m | |
|55 | 4078912619109602 | 485634974964502  | 45m | |


```
$ g++ -O3 -march=native -mavx2 -std=c++17 A000047.cpp ../utils/count_special_primes.cpp -lprimesieve -fopenmp -oA000047
$ time seq 20 40 | xargs -L1 ./A000047
...
| 40 | 146795295349       | 20601560352        | 1.65     |
A000047(40) = 146795295349

real	0m4.516s
user	0m14.454s
```

Python
```
$ time pypy A000047.py 31
	 52549929
A000047_final    (31) = 327770275

real	0m1.566s

$ time pypy A000047.py 40
	primes(82101) 2 ... [1049537, 1049549, 1049569]
	count_special_primes(2^40) = 20601560352
A000047_final    (40) = 146795295349  (60.04)

real	1m0.106s
```
