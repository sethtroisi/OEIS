# [A000049](https://oeis.org/A000049)

This sequence only had 36 terms (in 2022) extending it to 50 unblocked various OEIS self referential sequences.

See [A051070](https://oeis.org/A051070) and
[Self-referential sequences (wikipedia)](https://en.wikipedia.org/wiki/On-Line_Encyclopedia_of_Integer_Sequences#Self-referential_sequences)

## Results

| n  | a(n)          | enumerated    | queue time(s)  | hash time(s) |
|----|---------------|---------------|----------------|--------------|
| 25 | 3376376       | 7610774       | 0.45    secs   | |
| 26 | 6607207       | 15219614      | 0.95    secs   | |
| 27 | 12941838      | 30436779      | 1.96    secs   | |
| 28 | 25371086      | 60869794      | 4.07    secs   | |
| 29 | 49776945      | 121734504     | 8.39    secs   | |
| 30 | 97730393      | 243461570     | 17.37   secs   | |
| 31 | 192009517     | 486913005     | 34.05   secs   | |
| 32 | 377473965     | 973811315     | 71.07   secs   | |
| 33 | 742511438     | 1947601976    | 149.91  secs   | |
| 34 | 1461351029    | 3895174885    | 316.18  secs   | |
| 35 | 2877568839    | 7790308635    | 667.83  secs   | |
| 36 | 5668961811    | 15580558885   | 1405.92 secs   | |
| 37 | 11173141315   | 31161035403   | 2876.73 secs   | |
| 38 | 22030773337   | 62321953857   | 5939.09 secs   | |
| 39 | 43456681698   | 124643742959  | 12225.93 secs  | |
| 40 | 85752398532   | 249287251311  |                | 65/490 |
| 41 | 169273983901  | 498574171848  |                | 130/970 (3.826e9 iter/s) |
| 42 | 334256856592  | 997147875376  |                | 276/2124 |
| 43 | 660251251115  | 1994295089181 |                | 637/4940 |
| 44 | 1304578143057 | 3988589241254 |                | 1335/10377 (3.010e9 iter/s) |
| 45 | 2578437277523 | 7977177159689 |                | 3475/27095 (2.295e9 iter/s) |
| 46 | 5097564924796 | 15954352447918 |               | 13850/108094 |
| 47 | 10080525881679 | 31908702249586 |              | 42050/331300 |



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


| File/Method | Description |
|-------------|-------------|
| A000049_hash.cpp | Break into `num_class` congruence classes and shove in a large hash set |
  A000049_segmented_hash.cpp | Same as hash but does `num_passes`passes over each congruence class to minimize size of hash set |
| A000049_queue.cpp | Keep a priority queue of `(3*a^2 + 4*b^2, a, b)` at each step count and replace `(n, a, b)` with `(new_n, a, b+1)` |
| A000049_hash_with_queue.cpp | Combines hash and queue, breaking the problem into `num_class` classes and expanding each class with a queue` |


```
$ g++ -O3 A000049_queue.cpp && time ./a.out
$ g++ -fopenmp -O3 -std=c++17 -Wall -Werror A000049_hash.cpp && time ./a.out 34
$ g++ -O3 -fopenmp -Wall -Werror -std=c++17 A000049_hash_with_queue.cpp && time ./a.out 34
$ g++ -O3 -march=native -fopenmp -Wall -Werror -std=c++17 A000049_segmented_hash.cpp
```


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
