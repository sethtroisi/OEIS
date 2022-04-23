# [A000049](https://oeis.org/A000049)

This sequence currently (2022) only has 36 terms extending it to 50 unblocks various OEIS self referential sequences.

See [A051070](https://oeis.org/A051070) and
[Self-referential sequences (wikipedia)](https://en.wikipedia.org/wiki/On-Line_Encyclopedia_of_Integer_Sequences#Self-referential_sequences)

## Results
----------

| n | a(n) | iters | queue time(s)  | hash time(s) |
|---|------|-------|----------------|--------------|
| 26 | 6607207     | 15219613      |  0.89 secs | |
| 27 | 12941838    | 30436778      |  1.84 secs | |
| 28 | 25371086    | 60869793      |  3.79 secs | |
| 29 | 49776945    | 121734503     |  7.96 secs | |
| 30 | 97730393    | 243461569     | 16.62 secs | 6  |
| 31 | 192009517   | 486913005     | 34.05 secs | 12 |
| 32 | 377473965   | 973811315     | 71.07 secs | 21 |
| 33 | 742511438   | 1947601976    | 149.91 secs | 50 or 20/65 |
| 34 | 1461351029  | 3895174885    | 316.18 secs | 27/110 (elapse/total) |
| 35 | 2877568839  | 7790308635    | 667.83 secs | 50/270 |
| 36 | 5668961811  | 15580558885   | 1405.92 secs | 160/1160 |
| 37 | 11173141315 |    | 2876.73 secs |
| 38 | 22030773337 |    | 5939.09 secs |
| 39 | 43456681698 |    | 12225.93 secs |


```

| Method | Iterations / second (million) | Params |
|--------|-------------------------------|--------|
| Hash -> UnorderedMap       | 44-45 | 30-32, lots of classes |
| Hash -> ska::flat_hash_map |
| Hash -> Vector             | 20    | 30-31, lots of classes |
| Hash -> PriorityQueue      | 10-15 | 30-32, lots of classes |
| Queue                      | 11-16 | 28-36 |


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

1. [X] Hash approach can be trivially parallelized

Compared to the Queue approach this takes an older `O(log(C))` less time as hashset insert is `O(1)`.
In practice I believe they are about the same.

1. It's possible that shoving in a list and sorting would be even faster.
  * This is not true as for each unique item there's roughly 3x duplicates

## Ideas

* Can I split with a second modulo inside of each congruence class?
  * Probably not because of we're expanding `(x + i*base)`
* Can I skip all pairs where x and y share a factor?
* Can combine Queue Approach inside of Hash approach is hash table shows down too much
  * Tried this earlier with odd / even split and it didn't give any speed up
* All enumeration ideas (Hash / Queue) require enumerating and inserting on the order of 2^48 = 2.8e14 pairs

