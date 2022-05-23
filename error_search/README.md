# Searching for errors in OEIS

See

* http://oeis.org/wiki/Things_to_do_on_the_OEIS#Errors_to_correct
* https://davidbieber.com/snippets/2020-06-28-oeis-download/
* http://sequencedb.net/about.html

TODO

* A few more monotonic except 1 searches
* TODO search for old post about searching for all but one prime / automatically running code
* Consider using https://github.com/sidneycadot/oeis

Uses data from

* https://oeis.org/wiki/Welcome#Compressed_Versions


### To Inspect

* https://oeis.org/A145985 -> 753 seem wrong
* https://oeis.org/A177378 -> 67108861 seem wrong
* https://oeis.org/A195258 -> 7173 seem wrong (and probably related sequence https://oeis.org/A187825)
* https://oeis.org/A212605 -> 351439 seem wrong
* https://oeis.org/A228851 -> 1037917517 seem wrong (and probably related sequences)
* https://oeis.org/A248526 -> 53687090 seem wrong
* https://oeis.org/A279192 -> 692 seem wrong


### Found Manually

* https://oeis.org/A000022 -> Mathematica code disagrees at 50 (560915929775897218 vs 559792108243652284)
    * Jean-François Alcover's Mathematica program matches NJA Sloane's

* https://oeis.org/A228988 - Python code should be read(2, 10 ** a + 2), Mathematica code misses a digit too

* https://oeis.org/A003371 - Example is for https://oeis.org/A003367

* https://oeis.org/A210641 - a(59) is negative, update comment
* https://oeis.org/A000223 - a(59) is negative
    * Update sequence to give |P(A000092(n))|?
    * Update title to sequence gives nearest integer to P(A000092(n)) to match A000036
    * Could update to faster Mathematica
      ```
nmax = 5*10^4;
P[an_, n_] := an - Round[(4/3)*Pi*n^(3/2)];
record = 0;
an = 0;
A000223 =
 Reap[For[n = 0, n <= nmax, n++, an += SquaresR[3, n];
    If[(p = Abs[pn = P[an, n]]) > record, record = p; Print[pn];
     Sow[pn]]]][[2, 1]]
      ```


### To Fix

* https://oeis.org/A126788 -> Add b file and code
    * ```
      Primorial[n_] := Product[Prime[i], {i, 1, n}]
      A126788[n_] := (q = 1; While[!(PrimeQ[Prime[q]*Primorial[n] - 1] &&
            PrimeQ[Prime[q]*Primorial[n] + 1]), q++]; Return[Prime[q]])
      Table[A126788[n], {n, 1, 170}]
      {2, 2, 2, 2, 5, 11, 17, 11, 11, 37, 61, 23, 127, 53, 37, 1427, 23,
      11, 491, 11, 1429, 139, 83, 547, 233, 881, 149, 47, 313, 2, 463, 857,
      1367, 2269, 2221, 8171, 317, 619, 3727, 227, 173, 2251, 1091, 13,
      439, 277, 1597, 433, 173, 1303, 2927, 1367, 6679, 4903, 1979, 5641,
      1871, 11593, 3037, 1439, 5531, 17939, 823, 1381, 11677, 1861, 4591,
      9181, 293, 12391, 20113, 983, 701, 39791, 24691, 3607, 3331, 1087,
      51341, 1129, 28859, 4409, 18731, 14929, 28807, 1213, 3697, 10067,
      769, 4861, 3067, 34351, 16879, 17929, 521, 4507, 673, 46549, 2293,
      42131, 15031, 26597, 18089, 131, 4729, 42169, 10853, 10331, 577,
      11941, 919, 41813, 50387, 1051, 33311, 34961, 7723, 9011, 3739, 3671,
      53201, 461, 319439, 12497, 17551, 42557, 6883, 85607, 60289, 122207,
      30047, 23669, 10243, 45817, 6659, 43691, 2621, 4517, 35863, 114679,
      63977, 31033, 6529, 14627, 21481, 31193, 44497, 2, 119737, 34649,
      83059, 72077, 15749, 1291, 32141, 110291, 38377, 14563, 36563, 13697,
      22291, 97879, 55663, 31397, 93419, 103769, 971, 58727, 8317, 9343,
      16831, 156227, 36529, 49331, 260201, 5689, 55291, 17909, 62467, 1459,
      129263, 150961, 58211, 144611, 103, 42793, 38329, 36241, 1367, 44041,
      43319, 9473, 154417, 116293, 58757, 28513, 185291, 137353, 80713,
      48779}
      ```
* https://oeis.org/A172514 -> 2223344069957 is not prime should be 76375900241
    * Verified a(1 ... 25) not sure what to do about a(26)
    ```
an = {3, 7, 19, 97, 823, 3499, 2777, 6827, 2437, 21523, 300299, 446273, 339769, 1168523, 14117417, 29227421, 14160061, 78521987, 161187707, 1200085823, 2125209127, 1369430897, 56378083771, 26054006611};
ATestQ[n_, an_] := Block[{mult = n^(Floor[Log[n, an] + 2])}, !  AnyTrue[Table[l*mult + an*n + t, {l, 1, n - 1}, {t, 0, n - 1}], PrimeQ, 2]]
Table[SelectFirst[Prime /@ Range[1, PrimePi[an[[n - 1]]] + 10], ATestQ[n, #] &], {n, 2, 14}]
PrimeQ[2223344069957]
ATestQ[26, 2223344069957]
ATestQ[26, 76375900241]
    ```

* https://oeis.org/A209296 -> 670302913 should be 6703029313
    * Could add a(11) = 19161946762059695114295799175229653460447642068835611619471928021333660276190881810892144113
    * Mathematica `an = 1; Select[PrimeQ][Table[an += Times @@ Take[Select[Prime /@ Range[k + 20], ! Divisible[an, #] &], k], {k, 1, 20}]]`
* https://oeis.org/A246637 -> last term should be removed (or replaced with very large term)
* https://oeis.org/A257110 -> Prime sequence contained a composite (48315633)
    * Can add a bfile if so desired
* https://oeis.org/A283354 -> 71544139714543458911 should be 71544139714543459811
    * I think the next term should be 23, 838726639879696331, 115106274626071253862064364159, 3797, 167243045081386894442973727680539, 291781699730827438391998855229567883524768163044366987879

### Fixed

* https://oeis.org/A090208 -> 39088169 should be 102334155
    * `Total /@ Partition[ Select[Fibonacci /@ Range[0, 40], ! PrimeQ[#] &], 2, 1]`
* https://oeis.org/A114095 -> != should be == in Mathematica code
* https://oeis.org/A114831 -> 174392292 should be 1432401097
* https://oeis.org/A126788 -> 51 should be 61
* https://oeis.org/A136154 -> 16380 transposed error, should be 16830
* https://oeis.org/A138000 -> 1797439359 should be 1797439367
    * [ ] Could add 57518059831 if we want
* https://oeis.org/A138715 -> 3688 should be 36887
* https://oeis.org/A185656 -> 166778433667 should be 166778433637
    * Can also add 4054944029077, 58520914386689
* https://oeis.org/A237579 -> 2564940997072 should be 25649409970727
* https://oeis.org/A268465 -> 252077 is missing a final 9
* https://oeis.org/A307474 -> 17857 should be 178571
* https://oeis.org/A327914 -> 901 should be 907
* https://oeis.org/A340281 -> 590490 should be 5904901

Found Manually

* https://math.stackexchange.com/questions/2193384/i-think-i-found-an-error-in-a-oeis-sequence-what-is-the-proper-site-to-post-it

