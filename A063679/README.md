# [A063679](https://oeis.org/A063679)

Terms are 1-1 with a term in [A063680](https://oeis.org/A063680) and related to work on [A217259](../A217259).

`sieve.py` is a poor substitute for mtsieve.

After making a small
[changes](https://www.mersenneforum.org/showpost.php?p=616625&postcount=741)
to mtsieve I got it running and it's only 10-100x faster than sieve.py!

`./srsieve2 -fA -p 8 -P 4e8 -n 1000 -N 150000 -s "1*9^n-7"`

Ideally I'd use the OpenCL/GPU version but segfaults for unknown reasons.

## Results

I initially filtered N <= 300K to 34,436 terms with primes < 1e9 using `sieve.py`.
But my first pass of `sieve.py` didn't account for order(3, p). Then I discovered
`srsieve2` which is a much more mature tool.

Filtered N <= 400K to 23,168 terms with primes < 1e11 using srsieve2[^1]

Filtered 400K <= N <= 1M to 29,403 terms with primes < 1e13 using srsieve2[^1].
This took 10 core-days + 1 gpu-day.
srsieve's quoted removal rate was 334 seconds per factor but final delta was
closer to 600 seconds per factor.

[^1]: I had to remove some constraits and found some bugs/optimizations while doing so.
      See [#756 to #775](https://www.mersenneforum.org/showpost.php?p=617956&postcount=756)

```
# HAND MODIFY .abc file "1*9^$1-7" -> "(9^$1-7)/2"

# Starts pfgw64 per prime (OLD)
# tail +2 b9_n.abc | awk '{ print "(3^" $1 "-7)/2" }' | parallel -j 8 --lb 'eval pfgw64 -f0 -k -q"{}"'

# Split round-robin to X files
$ F=b9_huge_1e13.abc
$ tail +2 "$F" | split -d -n r/<X> - test_split. --filter='sh -c "{ head -n1 '$F'; cat; } > $FILE"'
# Call pfgw64 on each file
$ ls test_split.* | time taskset -c 0-9 parallel --lb -j 10 'pfgw64 -k -f1'

# Sort results (simple)
$ cat pfgw-prime.log pfgw.log | sort -u | sort -n -k1.4
# When results contain 9^exp
$ cat pfgw-prime.log pfgw.log | sed 's/(1\*/(/' | awk '/9\^/{match($0, /^\(9\^([0-9]*)-7\)\/2$/, t); print "(3^" 2*t[1] "-7)/2" ;}' | sort -u | sort -n -k1.4
```

Tested with parallel and pfgw64 to 3^500,000 over 48 hours and found a(14)-a(25)

```
(3^15398-7)/2
(3^15490-7)/2
(3^20408-7)/2
(3^39240-7)/2 <- Initially incorrectly reported, by me, as 34994
(3^41060-7)/2
(3^41842-7)/2
(3^58358-7)/2
(3^60346-7)/2
(3^82214-7)/2
(3^134972-7)/2
(3^194014-7)/2
(3^344204-7)/2
```


### Timing

Using 10 cores on a Ryzen 3900x

```
(9^200042-7)/2 is composite: RES64: [A6288AFAEDA86DC6] (77.2164s+0.4222s)
(9^220002-7)/2 is composite: RES64: [8CD532386BB8B894] (110.0737s+0.3480s)
```

