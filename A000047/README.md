# [A000047](https://oeis.org/A000047)

Quoting [A035251](https://oeis.org/A035251)

> A positive number n is representable in the form x^2 - 2y^2 iff every prime p == 3 or 5 (mod 8) dividing n occurs to an even power.

When I started this sequence only had 37 terms, extending it to 50 unblocks various OEIS self referential sequences.

See [A051070](https://oeis.org/A051070) and
[Self-referential sequences (wikipedia)](https://en.wikipedia.org/wiki/On-Line_Encyclopedia_of_Integer_Sequences#Self-referential_sequences)

I've adapted Lucy\_Hedgehog's [Count Primes](https://math.stackexchange.com/a/2283829/87805) to count the number of primes p == 3 or 5 (mod 8) under `n, n/2, n/3, ... ,n/sqrt(n), ... 5, 4, 3, 2, 1`  then I use inclusion-exclusion to

* I tried speeding up with `absl::flat_hash_map` and `ska::flat_hash_map` but it was 2x slower.
* One time in the past I tried converting `counts` to two arrays but the indexing was terribly hard.
  * **NEW:** I added a backing array so that I can iterate `V` and `counts[v]`, which gaves me a 2x speedup!


## Results
----------

|n  |count primes|A000047(n)|time|
|---|------------|----------|----|
|38 | 5433137518  | 37694785327 | 3s (c++) 23s (python) |
|39 | 10575969498 | 74373523152 | 5s (c++) 34s (python) |
|40 | 20601560352 | 146795295349 | 8s (c++) 62s (python) |
|41 | 40158311317 | 289837406442 | 14s (c++) 87s (python) |
|42 | 78330556749 | 572450196691 | 25s |
|43 | 152880952819 | 1130980093866 | 42s |
|44 | 298558223205 | 2235114244980 | 71s |
|45 | 583373483173 | 4418409686942 | 124s |
|46 | 1140499468654 | 8736714417829 | 210s (c++) 23m (python) |
|47 | 2230816679254 | 17279890334179 | 6.5m (c++) 33m (python) |
|48 | 4365594817954 | 34185319295647 | 11.4m (c++) <107m (python) |
|49 | 8547216668701 | 67645602070024 | 22.5m (c++) 132m (python)
|50 | 16741689949917 | 133886409875486 | 38m (c++) <262m (python) |
```
$ time seq 20 40 | xargs -I {} ./A000047 {}
...
A000047(40) = 146795295349

real	0m20.372s
```

C++
```
$ g++ -DNDEBUG -std=c++17 -O3 A000047.cpp -l primesieve -o A000047
$ time ./A000047 31
	count_special_primes(2^31) = 52549929
A000047(31) = 327770275

real	0m0.148s

$ time ./A000047 40
	count_special_primes(2^40) = 20601560352
A000047(40) = 153010765645

real	0m26.473s
```

Python
```
$ time pypy A000047.py 31
	 52549929
A000047_final    (31) = 327770275

real	0m1.732s

$ time pypy A000047.py 40
	 20601560352
A000047_final    (40) = 146795295349

real	1m32.377s
```
