import sys
sys.path.append('/usr/local/google/home/sethtroisi/Scripts/OEIS/')
import bfile
import random

import gmpy2
from gmpy2 import mpz
import MathLib

MAX_DIGITS = 1500

# Found
#   330 7 * 10^5880 + 3
# in 250 minutes = 15000s

seq = []

def LucasTest(a, d, b):
  assert b == 1
  primes = [2,5]
  for p in [3, 7]:
    if a % p == 0:
      primes.append(p)

  n = a * mpz(10) ** d + b
  for k in range(20):
    a = mpz(random.randint(2, 10000))
    if gmpy2.powmod(a, n - 1, n) != 1:
      return False
    for p in primes:
      if gmpy2.powmod(a, (n - 1) // p, n) == 1:
        break
    else:
      return True

  # probably composite
  return False

def addToSeq(a):
  t = str(a)
  prettyA = "{} * 10^{} + {}".format(t[0], len(t)-1, t[-1])
  seq.append(prettyA)

  print (len(seq), prettyA)

from collections import defaultdict
seen = defaultdict(int)

for p in range(1, 10):
  if gmpy2.is_prime(p):
    addToSeq(p)

leading_one = 1
for d in range(1, MAX_DIGITS):
  leading_one *= 10

  leading_digit = 0
  for a in range(1, 10):
    leading_digit += leading_one

    for b in [1]:
#    for b in range(1, 10, 2):
      if b == 5 or (a + b) % 3 == 0 or (a == b and d > 10):
        continue

      t = leading_digit + b

      if d > 10 and b == 1:
        is_prime = LucasTest(a, d, b)
#        assert lucas_is_prime == is_prime
      else:
        is_prime = gmpy2.is_prime(t)

      if is_prime:
        addToSeq(t)

#bfile.WriteListToFile("069675", seq)
