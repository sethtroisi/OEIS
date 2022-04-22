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
| 33 | 1280749500| 421577001 | 152.86 |

```


## Code

Density is 10% at 1M, 9.7% at 1B

64GB of RAM * 75% load factor / (int64 + 2 bytes overhead) = 4.8 billion
log2(4.8 / 0.09) = 36

So with giant brute force set can probably get to around existing terms.

---

Alternatively I can borrow from my taxicab strategy:

Keep a bunch of queues `(x, y)` with a fixed `x`.

Then use a priority queue to pull the smallest term `z = 2 x^2 + 4 y^2` and replace it in the queue with `(x, y+1)`


