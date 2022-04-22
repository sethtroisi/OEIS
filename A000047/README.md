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

Slow part is Lucy's Count Prime, I tried speeding up with `absl::flat_hash_map` but it was 2x slower.

One time in the past I tried converting `counts` to two arrays but the indexing was terribly hard.
