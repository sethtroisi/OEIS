# [A000049](https://oeis.org/A000049)

This sequence currently (2022) only has 36 terms extending it to 50 unblocks various OEIS self referential sequences.

See [A051070](https://oeis.org/A051070) and
[Self-referential sequences (wikipedia)](https://en.wikipedia.org/wiki/On-Line_Encyclopedia_of_Integer_Sequences#Self-referential_sequences)

## Results

| n  | a(n)         | iters         | queue time(s)  | hash time(s) |
|----|--------------|---------------|----------------|--------------|
| 25 | 3376376      | 7610774       | 0.45    secs   | |
| 26 | 6607207      | 15219614      | 0.95    secs   | |
| 27 | 12941838     | 30436779      | 1.96    secs   | |
| 28 | 25371086     | 60869794      | 4.07    secs   | |
| 29 | 49776945     | 121734504     | 8.39    secs   | |
| 30 | 97730393     | 243461570     | 17.37   secs   | 6, 19 |
| 31 | 192009517    | 486913005     | 34.05   secs   | 12, 38 |
| 32 | 377473965    | 973811315     | 71.07   secs   | 21, 80 |
| 33 | 742511438    | 1947601976    | 149.91  secs   | 50 or 20/65, 21/165, 30 |
| 34 | 1461351029   | 3895174885    | 316.18  secs   | 27/110 (elapse/total), 44/345, 63 |
| 35 | 2877568839   | 7790308635    | 667.83  secs   | 50/270, 94/730, 123 |
| 36 | 5668961811   | 15580558885   | 1405.92 secs   | 160/1160, 190/1525,  |
| 37 | 11173141315  | 31161035403   | 2876.73 secs   | 540/4200, 400/3168, 140/1140 |
| 38 | 22030773337  | 62321953857   | 5939.09 secs   | , 840/6660, 300/2500 |
| 39 | 43456681698  | 124643742959  | 12225.93 secs  | , 1770 / 14100, 700/5400 |
| 40 | 85752398532  | 249287251311  | ?              | , 4066 / 32404, 1800/12300 |
| 41 | 169273983901 | 498574171848  |                | , 10130 / 80550, 5400/42700 |


| Method | Iterations / second (million) | Params |
|--------|-------------------------------|--------|
| Hash -> `unordered_map`               | 44-45 | 30-32, lots of classes |
| Hash -> `ska::flat_hash_map`          | 40-50 | 30-34, lots of classes |
| Hash -> `vector` + sort + unique      | 20    | 30-31, lots of classes |
| Hash -> `priority_queue`              | 10-15 | 30-32, lots of classes |
| Queue -> `priority_queue`             | 11-16 | 28-36 |
| Queue -> `rollbear::prio_queue`       | 14-16 | 28-36 |
| SegmentedHash -> `ska::flat_hash_map` | 60+   | 33 (10007, 5 passes) |


```
$ g++ -O3 A000049_queue.cpp && time ./a.out
$ g++ -fopenmp -O3 -std=c++17 -Wall -Werror A000049_hash.cpp && time ./a.out 34
$ g++ -O3 -fopenmp -Wall -Werror -std=c++17 A000049_hash_with_queue.cpp && time ./a.out 34
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

1. [X] Hash approach can be trivially parallelized

Compared to the Queue approach this takes an older `O(log(C))` less time as hashset insert is `O(1)`.
In practice I believe they are about the same.

1. It's possible that shoving in a list and sorting would be even faster.
  * This is not true as for each unique item there's roughly 3x duplicates

## Ideas

* Can I split with a second modulo inside of each congruence class?
  * Probably not because of we're expanding `(x + i*base)`
* Can I skip all pairs where x and y share a factor?
* All enumeration ideas (Hash / Queue) require enumerating and inserting on the order of 2^48 = 2.8e14 pairs
* The Hash appoarch isn't doing anything productive with order of enumeration.
  Possibly I can iterate a small band of leading edge `(a, b)` for `b in [b_min, b_max]`
  This might take a little bit of overhead but would make the percent removed in each pass
  over the set much higher.

```
(a+c)^2 - a =
        2*a*c + c*c

Something's going on with
(112 | 2, 5)
(112 | 4, 4)
x gains 36
y loses 36

something about right triangles?
        one side is multiple of 3
        one side is multiple of 4
                36 / 3 = 12 => delta in a square progression (2 -> 4)
                36 / 4 = 9 => delta in a square progression (4 -> 5)

        choosing a random number divisible by 12
                72 / 3 = 24
                72 / 4 = 18
        when are 24 and 18 delta in square progressions?
                (gap must be even sized)
                24 = 2*a*c + c*c
                        c=2 => a=5
                                sure enought 7^2 - 5^2 = 24
                        c=4 => a=1
                                5^2 - 1^2 = 24
                18 = 2*a*c + c*c
                        c=2 => No
                        c=4 => No

A000050 / A000074
  * It's possible we can use prime signatures and subtract overcounting, See A057653
  * Also possible I can completely re-use A000047 code with special_primes being 4k+3
