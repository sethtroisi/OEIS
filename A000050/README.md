# [A000049](https://oeis.org/A000049)

This sequence only had 35 terms extending it to 50 unblocks various OEIS self referential sequences.

See [A051070](https://oeis.org/A051070) and
[Self-referential sequences (wikipedia)](https://en.wikipedia.org/wiki/On-Line_Encyclopedia_of_Integer_Sequences#Self-referential_sequences)

## Results

| n  | a(n)          | iters          | queue time(s)| hash time(s) |
|----|---------------|----------------|--------------|--------------|
| 25 | 6402706       | 13181716       | 0.8     secs | |
| 26 | 12534812      | 26360545       | 1.7     secs | |
| 27 | 24561934      | 52717056       | 3.6     secs | |
| 28 | 48168461      | 105428281      | 7.5     secs | |
| 29 | 94534626      | 210848467      | 15.5    secs | |
| 30 | 185661958     | 421685278      | 31.8    secs | |
| 31 | 364869032     | 843354269      | 65.9    secs | |
| 32 | 717484560     | 1686685460     | 138.0   secs | |
| 33 | 1411667114    | 3373338369     | 290.5   secs | |


| Method | Iterations / second (million) | Params |
|--------|-------------------------------|--------|
| Queue -> `priority_queue` | 11-16 | 28-36 |
| SegmentedHash -> `bitset` | 500+  | 37, bitset<4M> |

```
$ g++ -O3 A000050_queue.cpp && time ./a.out
$ g++ -O3 -march=native -fopenmp -Wall -Werror -std=c++17 A000050_segmented_hash.cpp
```

## Code

This code is borrowed near verbatim from [OEIS/A000049](../A000049).
