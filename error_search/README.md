# Searching for errors in OEIS

See

* https://math.stackexchange.com/questions/2193384/i-think-i-found-an-error-in-a-oeis-sequence-what-is-the-proper-site-to-post-it
* http://oeis.org/wiki/Things_to_do_on_the_OEIS#Errors_to_correct
* https://davidbieber.com/snippets/2020-06-28-oeis-download/

Inspect

* https://oeis.org/A145985 -> 753 seem wrong
* https://oeis.org/A177378 -> 67108861 seem wrong
* https://oeis.org/A195258 -> 7173 seem wrong (and probably related sequence https://oeis.org/A187825)
* https://oeis.org/A212605 -> 351439 seem wrong
* https://oeis.org/A228851 -> 1037917517 seem wrong (and probably related sequences)
* https://oeis.org/A248526 -> 53687090 seem wrong
* https://oeis.org/A279192 -> 692 seem wrong

Manual finds

* https://oeis.org/A000022 -> Mathematica code disagrees at 50 (560915929775897218 vs 559792108243652284)
  * Jean-François Alcover's Mathematica program matches NJA Sloane's

* https://oeis.org/A228988 - Python code should be read(2, 10 ** a + 2), Mathematica code misses a digit too

* https://oeis.org/A003371 - Example is for https://oeis.org/A003367


TODO

* https://oeis.org/A126788 -> 51 should be 61
    * ```
      Primorial[n_] := Product[Prime[i], {i, 1, n}]
      A126788[n_] := (q = 1; While[!(PrimeQ[Prime[q]*Primorial[n] - 1] &&
            PrimeQ[Prime[q]*Primorial[n] + 1]), q++]; Return[Prime[q]])
      Table[A126788[n], {n, 1, 1000}]
      '''
* https://oeis.org/A138000 -> 1797439359 should be 1797439367
* https://oeis.org/A237579 -> 2564940997072 should be 25649409970727
* https://oeis.org/A257110 - Prime sequence contained a composite (48315633)
* https://oeis.org/A327914 -> 901 should be 907

Fixed / Found

Opened a request for more extensions
