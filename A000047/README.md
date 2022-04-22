# [A000047](https://oeis.org/A000047)

Hoping to unblock the OEIS self referential sequence

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

$ $ time pypy A000047.py 40
	 20601560352
A000047_final    (40) = 146795295349

real	1m32.377s
```

* Slow part is Lucy's Count Prime, I tried speeding up with `absl::flat_hash_map` but it was 2x slower.
* One time in the past I tried converting `counts` to two arrays but the indexing was terribly hard.
  * **NEW:** I added a backing array so that I can iterate `V` and `counts[v]`, which gaves me a 2x speedup!


## Results
----------

|n  |count primes|A000047(n)|time|
|---|------------|----------|----|
|38 | 5433137518  | 37694785327 | 3s |
|39 | 10575969498 | 74373523152 | 5s |
|40 | 20601560352 | 146795295349 | 8s |
|41 | 40158311317 | 289837406442 | 14s |
|42 | 78330556749 | 572450196691 | 24s |
|43 | 152880952819 | 1130980093866 | 42s |
|44 | 298558223205 | 2235114244980 | 77s |
|45 | 583373483173 | 4418409686942 | 124s |
|46 | 1140499468654 | 8736714417829 | 210s (c++), 31m (python) |
|47 | 2230816679254 | 17279890334179 | 6.5m (c++) 46.5m (python) |
|48 | 4365594817954 | 34185319295647 | 11.4m (c++) 107m (python) |
|49 | 8547216668701 | 67645602070024 | 22.5m (c++)
|50 | 16741689949917 | 133886409875486 | 38m (c++) 262m (python) |
```
$ time seq 20 40 | xargs -I {} ./A000047 {}
...
A000047(40) = 146795295349

real	0m24.862s
```


