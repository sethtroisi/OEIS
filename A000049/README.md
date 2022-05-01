# [A000049](https://oeis.org/A000049)

This sequence only had 36 terms (in 2022) extending it to 50 unblocked various OEIS self referential sequences.

See [A051070](https://oeis.org/A051070) and
[Self-referential sequences (wikipedia)](https://en.wikipedia.org/wiki/On-Line_Encyclopedia_of_Integer_Sequences#Self-referential_sequences)

At somepoint I had at least one bug in my code (around handling (0,0) and the nodes it creates).
This lead to `count_pairs.py` which should be used to verify that the correct number of pairs are
enumerated. This check could be automated but isn't at this time.

No root cause was ever found and the issue was an overcounting in a(n) but an undercounting in
enumerated so it's a tiny bit scary! Sadly it only seemed to appear above a(45) so it's not easy
to isolate.

## Results

| 40 | 85752398532 | 249287251310  | 27673.14 secs | size: 605395, iters/s: 9.0 million

| n  | a(n)           | enumerated      | queue time(s)  | hash time(s) |
|----|----------------|-----------------|----------------|--------------|
| 25 | 3376376        | 7610773         | 0.45           | |
| 26 | 6607207        | 15219613        | 0.95           | |
| 27 | 12941838       | 30436778        | 1.96           | |
| 28 | 25371086       | 60869793        | 4.07           | |
| 29 | 49776945       | 121734503       | 8.39           | |
| 30 | 97730393       | 243461569       | 17             | |
| 31 | 192009517      | 486913005       | 34             | |
| 32 | 377473965      | 973811315       | 71             | |
| 33 | 742511438      | 1947601976      | 150            | |
| 34 | 1461351029     | 3895174885      | 316            | |
| 35 | 2877568839     | 7790308635      | 668            | |
| 36 | 5668961811     | 15580558885     | 1406           | |
| 37 | 11173141315    | 31161035403     | 2877           | |
| 38 | 22030773337    | 62321953856     | 5939           | 15/119 |
| 39 | 43456681698    | 124643742958    | 12225          | 29/229 |
| 40 | 85752398532    | 249287251310    | 27673          | 60/480 |
| 41 | 169273983901   | 498574171847    | 64008          | 130/970 (3.826e9 iter/s) |
| 42 | 334256856592   | 997147875375    |                | 270/2143 (3.686e9 iter/s) |
| 43 | 660251251115   | 1994295089542   |                | 610/4812 (3.266e9 iter/s) |
| 44 | 1304578143057  | 3988589241765   |                | 1335/10377 (3.010e9 iter/s) |
| 45 | 2578437277523  | 7977177160412   |                | 3540/28200 (2.254e9 iter/s) |
| 46 | 5097564924599  | 15954352449461  |                | 10304/82020 (1.548e9 iter/s) |
| 47 | 10080525879900 | 31908702253605  |                | 30197/240400 (1.056e9 iter/s) |


| File/Method | Description |
|-------------|-------------|
| A000049_hash.cpp | Break into `num_class` congruence classes and shove in a large hash set |
  A000049_segmented_hash.cpp | Same as hash but does `num_passes`passes over each congruence class to minimize size of hash set |
| A000049_queue.cpp | Keep a priority queue of `(3*a^2 + 4*b^2, a, b)` at each step count and replace `(n, a, b)` with `(new_n, a, b+1)` |
| A000049_hash_with_queue.cpp | Combines hash and queue, breaking the problem into `num_class` classes and expanding each class with a queue` |


```
# Fastest and best
$ g++ -O3 -march=native -fopenmp -Wall -Werror -std=c++17 A000049_segmented_hash.cpp
# Can help with cache contention
$ export OMP_PROC_BIND=spread
# Possibly 10% faster (only <= 42?)
$ clang++ -fopenmp=libiomp5 ...

# Double check / extra methods
$ g++ -O3 A000049_queue.cpp && time ./a.out
$ g++ -O3 -fopenmp A000049_hash.cpp && time ./a.out 34
$ g++ -O3 -fopenmp -std=c++17 A000049_hash_with_queue.cpp && time ./a.out 34
```


### Timing / Performance

Google Cloud Platform

| Platform                      | lscpu                      | vCPU | speed                  |
| e2-24core                     | Intel Xeon @ 2.2 (55MB L3) | 12/24| 40,16MB -> 87s, 2853   |
| c2-8                          | Intel Xeon @ 3.1 (24MB L3) | 4/8  | 40,16MB ->      1200   |
| c2-8                          | Intel Xeon @ 3.1 (24MB L3) | 4/8  | 40,12MB ->      1200   |
| c2d-standard-8 ("AMD Milan")  | AMD EPYC 7B13  (32MB L3)   | 4/8  | 40,32MB -> 134s, 1850  |
| c2d-standard-8 ("AMD Milan")  | AMD EPYC 7B13  (32MB L3)   | 4/8  | 40,16MB -> 128s, 1941  |
| c2d-standard-8 ("AMD Milan")  | AMD EPYC 7B13  (32MB L3)   | 4/8  | 41,32MB -> 268s, 1850  |

| n2d-standard-8 ("AMD Rome")   | AMD EPYC 7B12  (16MB L3)   | 8    | 40,32MB -> 197s, 1264  |
| n2d-standard-8 ("AMD Rome")   | AMD EPYC 7B12  (16MB L3)   | 8    | 40,16MB -> 165s, 1510  |
| t2d-standard-8 ("AMD Milan")  | AMD EPYC 7B13a (32MB L3)   | 8    | 40,32MB -> 94s, 2655   |
| t2d-standard-8 ("AMD Milan")  | AMD EPYC 7B13a (32MB L3)   | 8    | 40,16MB -> 77s, 3235   |
| t2d-standard-8 ("AMD Milan")  | AMD EPYC 7B13a (32MB L3)   | 8    | 41,32MB -> 166s, 3000  |
| t2d-standard-8 ("AMD Milan")  | AMD EPYC 7B13a (32MB L3)   | 8    | 41,16MB -> 146s, 3401  |

| t2d-standard-16 ("AMD Milan") | AMD EPYC 7B13a (32MB L3)   | 16   | 41,16MB -> 84s, 5935   |
| t2d-standard-16 ("AMD Milan") | AMD EPYC 7B13a (32MB L3)   | 16   | 41,32MB -> 92s, 5400   |

| t2d-standard-48 ("AMD Milan") | AMD EPYC 7B13a (32MB L3)   | 48   | 40,32MB -> 15s, 16350  |
| t2d-standard-48 ("AMD Milan") | AMD EPYC 7B13a (32MB L3)   | 48   | 41,32MB -> 31s, 15850  |


| Method | Iterations / second (million) | Params |
|--------|-------------------------------|--------|
| Hash -> `unordered_map`               | 44-45 | 30-32, lots of classes |
| Hash -> `ska::flat_hash_map`          | 40-50 | 30-34, lots of classes |
| Hash -> `vector` + sort + unique      | 20    | 30-31, lots of classes |
| Hash -> `priority_queue`              | 10-15 | 30-32, lots of classes |
| Queue -> `priority_queue`             | 11-16 | 28-36 |
| Queue -> `rollbear::prio_queue`       | 14-16 | 28-36 |
| SegmentedHash -> `ska::flat_hash_map` | 60+   | 33 (10007, 5 passes) |
| SegmentedHash -> `bitset`             | 500+  | 37, bitset<32M> |

## Queue approach

I reused an approach I learned at Rose-Hulman for enumerating taxicab numbers

A priority queue is filled with tuples `(3x^2 + 4y^2, x, y)` each with `x` fixed.

At each step pull the smallest item and replace it in the queue with `(x, y+1)`

This theoretically takes `O(log2(C))` time for each iteration with `C` slowly growing from 10 to ~200,000.

`log2(200,000) = 18`, Yet the code only manages 15 million iterations / second. So this approach has it's limitations.

* Can change to `4x^2 + 3y^2` which slightly reduces number of items in queue -> 1-3% speedup

* Technically it's possible to parallelize, would have each look at a range of numbers
  * for each `(x, y)` can easily find first `y` such that `n > start`

## Hash approach

1. Enumerate all pairs (possibly in 2 or 4 congruence classes), adding them to a hashmap.
1. Count number of items in the hashmap.

* [X] Hash approach can be trivially parallelized
* [X] Skip some initial passes by Radix sorting items (into which pass they first belong in)
* [X] (tested slower) Write `(n - pass_min) >> shift` to array in1st loop. In 2nd loop set `bitset.set(B[i])`

In theory hash set is `O(1)` but `unordered_set` and the rest all have a large constant
that makes `O(log(sqrt(N)))` Queue approach take less time in practice.

But then I remember the perfect Hash implementation is actually `bitset<>`.
In testing `bitset` is 20-50 times faster!

## Ideas

* Can I split with a second modulo inside of each congruence class?
  * Probably not because of we're expanding `(x + i * base)`
* Can I skip all pairs where x and y share a factor?
* Can I do something like odd / even inside of a congruence class to make bitset
  represent 2x more?

* Can I skip many passes after being included in a pass? (No)

With `num\_classes ~ 2^13, bitset<2 ^ 25 = 32MB>` and `n = 2^50, x = 0, y = 2^24`

```
y_delta = eight_base * y + four_base_squared
bitset represents 2^25 * 2^13

y_delta = 8 * 2^13 * 2^24 + 4 * (2^13)^2 = 8 * 2^13 * 2^24 = 2^40
y_delta / bitset ~ 8 * sqrt(N) / bitset = 8 * 2^24 / 2^25 = 4
```

If `sqrt(N) > bitset.size()` this would make more sense.

* What happens when we increase `num\_classes` by 2x?
  * 4x more setup work (in `build\_congruences`)
  * bitset represents 2x more -> half as many `num\_passes`
  * Seems to be mostly a wash?
    * More setup work
    * Few X expanded at a time (but readahead is probably great)
