# [A000049](https://oeis.org/A000049)

This sequence currently (2022) only has 36 terms extending it to 50 unblocks various OEIS self referential sequences.

See [A051070](https://oeis.org/A051070) and
[Self-referential sequences (wikipedia)](https://en.wikipedia.org/wiki/On-Line_Encyclopedia_of_Integer_Sequences#Self-referential_sequences)


## Results
----------

|n  | a(n) | iters | time |
|---|------|-------|------|
| 26 | 6607207 |    4022309 | 0.95 |
| 27 | 12941838|    7996869 | 1.98 |
| 28 | 25371086|   15894288 | 4.11 |
| 29 | 49776945|   31581134 | 8.50 |
| 30 | 97730393|   62738514 | 17.45 |
| 31 | 192009517|  124616045 | 35.78 |
| 32 | 377473965|  247490677 | 74.27 |
| 33 | 742511438 |  491483603 | 151.65 |
| 34 | 1461351029 |  975970220 |318.25 |
| 35 | 2877568839 | 1937994971 | 664.40 |
| 36 |  5668961811 | 3848260862 | 1386.74 |
| 37 |  11173141315 | 7641513003 | 2876.73 |
| 38 |  22030773337 | 15174081879 | 5939.09 |


```


## Queue approach

I reused an approach I learned at Rose-Hulman for enumerating taxicab numbers

A priority queue is filled with tuples `(3x^2 + 4y^3, x, y)` each with `x` fixed.

At each step pull the smallest item and replace it in the queue with `(x, y+1)`

This theoretically takes `O(log2(C))` time for each iteration with `C` slowly growing from 10 to ~200,000.

`log2(200,000) = 18`, Yet the code only manages 2.2 million iterations / second. So this appoarch has it's limitations.

## Hash approach

1. Enumerate all pairs (possibly in 2 or 4 congruence classes), adding them to a hashmap.
1. Count number of items in the hashmap.

Compared to the Queue approach this takes an older `O(log(C))` less time as hashset insert is `O(1)`.
In practice I believe they are about the same.

It's possible that shoving in a list and sorting would be even faster.

`64GB * 75% load factor / (64bits + 4byte overhead) = 2.5e9`.
The population has <10% density.
Would "only" need to break into `2^50 / 10 / 2.5e9` = 45035 classes.
2^20 classes would probably mops this up nicely.

## Ideas

1. Can split to even `(2*i, j)` and odd pairs `(2*i+1, j)` but that doesn't make the queue solution any faster.
  * From https://oeis.org/A001235/a001235_1.txt it's possible I can split into arbitary many classes
    * Pro: Makes memory demand manageable
    * Pro: Makes each hash insert can be faster.
    * Con: Still have to enumerate every item

1. Can I skip all pairs where x and y share a factor?

1. All enumeration ideas still require enumerating something on the order of 2^45 = 10^13.5 pairs.

