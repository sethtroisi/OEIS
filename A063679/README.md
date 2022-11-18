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

Tested `<= 3^500,000` with parallel and pfgw64 over 48 hours and found a(14)-a(25)

Tested `<= 3^1,000,000` with parallel (over 10 cores) over 10 days and found a(26)-a(27)

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
(3^587712-7)/2
(3^778070-7)/2
```


### Timing

Using 10 cores on a Ryzen 3900x

```
(9^200042-7)/2 is composite: RES64: [A6288AFAEDA86DC6] (77.2164s+0.4222s)
(9^220002-7)/2 is composite: RES64: [8CD532386BB8B894] (110.0737s+0.3480s)
(9^240019-7)/2 is composite: RES64: [145D3157E8D6FA37] (120.8113s+0.4579s)
(9^250012-7)/2 is composite: RES64: [6E1DBB488BED25C6] (122.1138s+0.4418s)
(9^260006-7)/2 is composite: RES64: [956C550FAE10F794] (130.3365s+0.4787s)
(9^280000-7)/2 is composite: RES64: [A306107F9960A6E9] (174.9219s+0.5521s)
(9^293856-7)/2 is 3-PRP! (193.7032s+0.6723s)
(9^300023-7)/2 is composite: RES64: [138BC72A8318E2F6] (185.5669s+0.6347s)
(9^320000-7)/2 is composite: RES64: [212FB5ABE6C9E478] (206.6544s+0.7203s)
(9^340010-7)/2 is composite: RES64: [7615C5C9F82784C2] (255.2589s+0.8084s)
(9^350008-7)/2 is composite: RES64: [041F8FABDEEB667C] (261.5549s+0.8535s)
(9^360007-7)/2 is composite: RES64: [C6EF464AB3ED4B33] (275.5616s+0.9180s)
(9^380030-7)/2 is composite: RES64: [629594DA01123884] (318.3652s+1.0155s)
(9^389035-7)/2 is 3-PRP! (325.2043s+1.0612s)
(9^400032-7)/2 is composite: RES64: [5ACB1BA039050C34] (332.5460s+1.1244s)
(9^450004-7)/2 is composite: RES64: [2C9EAF2D6BE84595] (449.3393s+1.4077s)
(9^499975-7)/2 is composite: RES64: [4814F4562B4F8ABC] (518.8192s+1.7371s)
```

Timing for intervals is fairly well estimated by
```python
>>> with open("b9_n_1e13.abc") as f: D = list(map(int, f.readlines()[1:]))
>>> R = range(200000, 240000); sum([4.1 ** math.log2(e / 200000) for e in D if e in R]) * 85 / (10 * 3600)
230 [hours]
```
