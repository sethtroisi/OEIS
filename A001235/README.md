# [A001235](https://oeis.org/A001235)

Fast enumeration of Taxi-cab numbers: sums of 2 cubse in more than 1 way

## Taxicab with Priority Queue
`taxicab_queue` uses a min `priority_queue` which keeps one item for each a value in a tuple `(a^3+b^3, a, b)`

At each step the smallest `(a^3 + b^3, a, b)` is replaced by `(a^3 + (b+1)^3, a, b+1)`.
Occasionally a new `(a^3, a, 0)` entry needs to be inserted.

This technique requires storing a few thousand tuples at a time which is quite fast in practice.
The throw-everything-in-a-hashset `O(1)` method is theoretically faster,
but has a much large constant for the ranges of interest.

## Taxicab with Hash Set

Rough code is done, see [A000049](../A000049/README.md) for two implemented ideas.

1. Break `a^3 + b^3` into congruence classes (which have no overlap)
1. Enumerate in several passes so fewer removal passes are needed over the HashSet

## Results

```
$ g++ -O3 taxicab_queue.cpp -o taxicab_queue
$ ./taxicab_queue
        1th           1729  (b:    1 size: 3  time: 0.0  iters: 74)
        2th           4104  (b:    2 size: 4  time: 0.0  iters: 128)
        3th          13832  (b:    2 size: 5  time: 0.0  iters: 275)
        4th          20683  (b:   19 size: 6  time: 0.0  iters: 359)
        5th          32832  (b:    4 size: 7  time: 0.0  iters: 482)
       10th          65728  (b:   12 size: 8  time: 0.0  iters: 758)
      100th        4673088  (b:   25 size: 35  time: 0.0  iters: 12501)
      200th       16387189  (b:    5 size: 53  time: 0.0  iters: 28711)
      500th      106243219  (b:  363 size: 98  time: 0.0  iters: 99496)
     1000th      416318561  (b:  105 size: 154  time: 0.0  iters: 246913)
     1100th      503842304  (b:  280 size: 164  time: 0.0  iters: 280363)
     1200th      599881464  (b:   93 size: 174  time: 0.0  iters: 314914)
     1500th      944972288  (b:  336 size: 203  time: 0.0  iters: 426204)
     2000th     1671816384  (b:  940 size: 244  time: 0.0  iters: 623168)
     4000th     7103146447  (b:  759 size: 397  time: 0.1  iters: 1633761)
     6000th    16541054656  (b: 1292 size: 525  time: 0.1  iters: 2869469)
     8000th    30225888875  (b: 1435 size: 643  time: 0.2  iters: 4288332)
    10000th    48538460952  (b: 1747 size: 752  time: 0.3  iters: 5880086)
    12000th    72125011153  (b: 1974 size: 859  time: 0.4  iters: 7656386)
    14000th   100547229384  (b: 3202 size: 960  time: 0.5  iters: 9554130)
    16000th   133882383096  (b: 3470 size: 1055  time: 0.6  iters: 11563099)
    18000th   172921387464  (b: 1122 size: 1150  time: 0.7  iters: 13713482)
    20000th   216264806875  (b:  930 size: 1238  time: 0.8  iters: 15918176)
    30000th   520890296211  (b:  195 size: 1660  time: 1.5  iters: 28599719)
    40000th   976889700163  (b: 6894 size: 2047  time: 2.3  iters: 43491829)
    50000th  1601017842872  (b: 2464 size: 2413  time: 3.2  iters: 60454343)
    60000th  2389857538048  (b: 10136 size: 2758  time: 4.3  iters: 78958931)
    70000th  3361947334704  (b: 7470 size: 3090  time: 5.5  iters: 99130503)
    80000th  4514221398024  (b: 6008 size: 3410  time: 6.8  iters: 120650540)
    90000th  5858794561158  (b: 9121 size: 3719  time: 8.2  iters: 143551914)
   100000th  7420456017000  (b: 14320 size: 4023  time: 9.7  iters: 168042671)
   110000th  9172388012625  (b: 3505 size: 4318  time: 11.3  iters: 193546111)
   120000th 11149741127168  (b:  732 size: 4609  time: 13.0  iters: 220446894)
   200000th 35059220195419  (b: 12622 size: 6752  time: 29.1  iters: 473132787)
   300000th 87566615553635  (b: 32858 size: 9161  time: 55.9  iters: 870963723)
   400000th 168270869224640  (b: 30208 size: 11390  time: 89.2  iters: 1346200578)
   500000th 280196391751731  (b: 50483 size: 13499  time: 128.1  iters: 1891223793)
   600000th 425167544990279  (b: 23191 size: 15512  time: 172.3  iters: 2497312491)
   700000th 604965422789784  (b: 55686 size: 17448  time: 221.6  iters: 3159259277)
   800000th 822617684187201  (b: 55908 size: 19330  time: 275.3  iters: 3877595955)
   900000th 1079224830750033  (b: 78705 size: 21161  time: 333.6  iters: 4646981439)
  1000000th 1376426637022528  (b: 9268 size: 22948  time: 396.4  iters: 5465095540)
```
