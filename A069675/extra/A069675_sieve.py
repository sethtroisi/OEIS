import sys
sys.path.append('/usr/local/google/home/sethtroisi/Scripts/OEIS/')
import bfile

import gmpy2
import time

import MathLib

from collections import defaultdict


T0 = time.time()

MAX_DIGITS = 2000
SIEVE_MAX = 5000

to_test = [[[True, True, False, True, True] for a in range(10)] for d in range(MAX_DIGITS+1)]

# Filter the 6 silly patterns
for a in (3, 6, 9):
  for b in (3, 9):
    for test_d in range(MAX_DIGITS + 1):
      assert (a * pow(10, test_d, 3) + b) % 3 == 0
      to_test[test_d][a][b // 2] = False

# Sieve out "small" prime factors and mark those numbers not to test.
small_primes = MathLib.sieveOfErat(SIEVE_MAX)
factors = MathLib.sieveOfFactors(SIEVE_MAX)

print ("\t PrimePi({}) = {}".format(SIEVE_MAX, len(small_primes)))

for p in small_primes:
  if p in (2, 5):
    continue

  divisible_mods = defaultdict(list)

  for a in range(1, 10):
    if a % p == 0:
      continue

    modular_inverse = MathLib.modularInverse(a, p)
    for b in (1,3,5,9):
      t = (p - b) * modular_inverse % p
      divisible_mods[t].append((a, b))

  order = 1
  power_ten = 10
  while power_ten != 1:
    order += 1
    power_ten = (10 * power_ten) % p

  assert pow(10, order, p) == 1

  power_ten = 1
  for d in range(max(order, MAX_DIGITS) + 1):
    if power_ten in divisible_mods:
      for a, b in divisible_mods[power_ten]:
        for test_d in range(d, MAX_DIGITS + 1, order):
          assert (a * pow(10, test_d, p) + b) % p == 0, "{} * 10 ^ {} + {} % {} != 0".format(a, d, b, p)
          to_test[test_d][a][b // 2] = False

    power_ten = (10 * power_ten) % p

filtered = sum(not a_v for d_v in to_test for b_v in d_v for a_v in b_v)
filtered -= 1 * 10 * (MAX_DIGITS + 1)
total_filterable = 4 * 10 * (MAX_DIGITS + 1)

print ()
print ("\t filtered {} of {} = {:.3f}".format(filtered, total_filterable, filtered / total_filterable))

T1 = time.time()

# Verify
# '''
print ()
for d in range(1, MAX_DIGITS + 1):
  for a in range(1, 10):
    for b in (1, 3, 5, 9):
      if to_test[d][a][b // 2]:
        t = a * 10 ** d + b
        for p in small_primes:
          if t % p == 0:
            print ("what what?", p, "\t", d, a, b)
# '''

T2 = time.time()

print ("\t filter: {:.1f}s (verify: {:.1f}s)".format(T1 - T0, T2 - T1))


seq = []

def addToSeq(a):
  seq.append(a)

  t = str(a)
  prettyA = "{} * 10^{} + {}".format(t[0], len(t)-1, b)
  if len(seq) % 13 == 0:
    print (len(seq), prettyA)

for p in range(1, 10):
  if gmpy2.is_prime(p):
    addToSeq(p)

leading_one = 1
for d in range(1, MAX_DIGITS):
  leading_one *= 10

  # (a,b) creates 9*4 = 36 pairs
  #   => turn pair into a module (of leading_one) that will be zero
  #     ex.
  #       p = 137, a = 2, b = 5
  #       2 * m + 5 = 0 (mod 137)
  #       2 * m = 132 (mod 137)
  #         (modular inverse of 2 mod 137 = 69 => 2 * 69 = 1 mod 137)
  #       69 * 2 * m = 69 * 132 (mod 137)
  #       m = 66
  #   => if leading_one in divisble_mods
  #   => mark d + k * p, a, b as not prime

  leading_digit = 0
  for a in range(1, 10):
    leading_digit += leading_one

    for b in (1,3,5,9):
      if to_test[d][a][b // 2]:
        t = leading_digit + b

        if gmpy2.is_prime(t):
          addToSeq(t)

bfile.WriteListToFile("069675", seq)

T3 = time.time()

print ("\t main: {:.1f}s".format(T3 - T2))
print ("\t total: {:.1f}s".format(T3 - T0))
