import gmpy2
import itertools
import subprocess
import math
import time


from factordb.factordb import FactorDB

from collections import Counter, defaultdict

# Also see A056938

def product(factors):
  temp = 1
  for factor in factors:
    temp *= factor
  return temp


def factor_large(n, b1=10**6):
  args = ["ecm", "-q", "-c", "10", str(b1)]
  print ("\t\t", " ".join(args))
  result = subprocess.run(
      args,
      input=str(n).encode(),
      stdout=subprocess.PIPE)

  if result.returncode == 8:
    # Need to rerun with smaller b1
    print("\t\tfound self ({} with b1={})".format(n, b1))
    return factor_large(n, b1= max(1000, b1 // 100))

  return list(map(int, result.stdout.strip().split()))


def attempt_factorization(s, know_factors):
  t = s
  factors = []

  for factor in know_factors:
    assert gmpy2.is_prime(factor)
    t //= factor
    factors.append(factor)

  #if False:
  if t >= 1e60:
    # Check factorDB (probably already been done)
    time.sleep(0.2)
    factordb = FactorDB(s)
    factordb.connect()
    factordb_factors = factordb.get_factor_list()
    if factordb_factors and factordb_factors != [s]:
      print ("\t\tfactordb:", factordb.get_status(), factordb_factors)
      for factor in factordb_factors:
        if gmpy2.is_prime(factor):
          t //= factor
          factors.append(factor)

  if t > 1 and t < 1e40:
    print ("\t\tpre:", factors, t)

    while t > 1:
      if gmpy2.is_prime(t):
        factors.append(int(t))
        t //= t
        break

      ecm_factors = sorted(factor_large(t))
      print ("\t\tecm:", ecm_factors)
      for factor in ecm_factors:
        if gmpy2.is_prime(factor):
          t //= factor
          factors.append(int(factor))

  return t, factors


def load_from_file():
  home_primes = defaultdict(list)

  n = None
  s = None
  with open("home_primes.txt") as f:
    # each line is "<base> <start> <step> <status>: <factor> <factor> ..."
    for line in f.readlines():
      pre, post = line.strip().split(":")
      *pre, status = pre.split()
      base, start, step, = map(int, pre)

      if start != n:
        n = start
        s = n

      factors = list(map(int, post.split()))
      if status not in ("FF", "P"):
        print (line)
      else:
        home_primes[(base, start, step)] = factors
        assert product(factors) == s, (start, step, s, factors)
        s = int("".join(map(str, factors)))


  all_primes = set()
  for factors in home_primes.values():
    for p in factors:
      if gmpy2.is_prime(p):
        all_primes.add(p)

  print ("Found", len(all_primes), "primes")
  return home_primes, all_primes


# For use with kernprof -v --line-by-line simple.py
#@profile
def run():
  home_primes, all_primes = load_from_file()

  required_steps = Counter()
  max_prime = Counter()
  try:
    for n in range(2, 100):
      print (n)
      t = n
      for step in itertools.count(0):
        if gmpy2.is_prime(t):
          break

        s = t

        original = home_primes[(10, n, step)]
        t, factors = attempt_factorization(s, original)

        if t > 1:
          print ("Breaking, failed to factor C{}: {}".format(len(str(t)), t))
          break

        factors.sort()
        assert product(factors) == s, (s, t, factors)

        new = int("".join(map(str, factors)))
        t = new

        if step >= 0 or gmpy2.is_prime(s):
          if new < 1e40:
            print ("\t", step, new, "from", s, factors)
          else:
            print ("\t", step, new)
            print ("\t\tfrom", factors)

        if factors != original:
          home_primes[(10, n, step)] = factors
          print ("\t\tnew factor", factors)

        # Size of 2nd largest prime
        max_prime[len(str(factors[-2:][0]))] += 1

      if gmpy2.is_prime(t):
        home_primes[(10, n, step)] = [t]
        required_steps[step] += 1
      else:
        required_steps[1000] += 1
        print ("\t {} Gave({}th time) up on step {}".format(n, required_steps[1000], step))

  except KeyboardInterrupt:
    print("Stopping from ^C")

  with open("home_primes.txt", "w") as f:
    for base, start, step in sorted(home_primes.keys()):
      factors = home_primes[(base, start, step)]
      if not factors:
        continue

      if all(gmpy2.is_prime(f) for f in factors):
        if len(factors) == 1:
          status = "P"
        else:
          status = "FF"
      else:
        status = "CF"

      f.write("{} {} {} {}: {}\n".format(
          base, start, step, status, " ".join(map(str, factors))))

  print ()
  for s, c in sorted(required_steps.most_common()):
    print ("{} with {} steps".format(c, s))
  print ("1000 steps = stopped early")

  print ()
  for p, c in sorted(max_prime.items()):
    print ("2nd largest primes had {} digits in {} steps".format(p, c))


run()
