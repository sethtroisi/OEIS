import gmpy2
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


def trial_factor(t, primes, original):
  t_sqrt = math.sqrt(t)

  factors = []
  for p in primes:
    if p > t_sqrt:
      break

    while t > 0 and t % p == 0:
      factors.append(p)
      t //= p
      t_sqrt = math.sqrt(t)

      if t in original or gmpy2.is_prime(t):
        factors.append(t)
        t //= t
        t_sqrt = 1

  return t, factors


def attempt_factorization(s, primes, original):
  # This is the slowest part of rerunning
  t, factors = trial_factor(s, primes, original)

  if t >= 1e40:
    # Check factorDB
    time.sleep(0.2)
    factordb = FactorDB(s)
    factordb.connect()
    factordb_factors = factordb.get_factor_list()
    if factordb_factors and factordb_factors != [s]:
      print ("\t\tfactordb:", factordb.get_status(), factordb_factors)
      if factordb.get_status() == "FF":
        assert set(factors) < set(factordb_factors)
        return 1, factordb_factors

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


# For use with kernprof -v --line-by-line simple.py
#@profile
def run():
  with open("factors") as f:
    extra_primes = sorted(map(int, f.readlines()))

  home_primes = defaultdict(list)
  with open("home_primes.txt") as f:
    # each line is "<base> <start> <step> <status>: <factor> <factor> ..."
    for line in f.readlines():
      pre, post = line.strip().split(":")
      *pre, status = pre.split()
      base, start, step, = map(int, pre)

      factors = list(map(int, post.split()))
      home_primes[(base, start, step)] = factors

      if status not in ("FF", "P"):
        print (line)

  # copy to see if we need to write new_factors
  original = set(extra_primes)

  print ("Found", len(extra_primes), "extra primes")
  for p in extra_primes:
    assert gmpy2.is_prime(p), p
  print ("Verified")

  primes = [2]
  while primes[-1] < 10 ** 4:
    primes.append(int(gmpy2.next_prime(primes[-1])))
  primes = sorted(set(primes) | set(extra_primes))

  required_steps = Counter()
  max_prime = Counter()
  try:
    for n in range(3500, 10000):
      print (n)
      step = 1
      t = n
      while not gmpy2.is_prime(t):
        s = t

        t, factors = attempt_factorization(t, primes, original)

        if t > 1:
          print ("Breaking, failed to factor C{}: {}".format(len(str(t)), t))
          break

        factors.sort()
        assert product(factors) == s, (s, t, factors)

        new = int("".join(map(str, factors)))
        t = new

        if step > 40 or gmpy2.is_prime(t):
          if new < 1e40:
            print ("\t", step, new, "from", s, factors)
          else:
            print ("\t", step, new)
            print ("\t\tfrom", factors)

        home_primes[(10, n, step-1)] = factors

        extra_primes += [f for f in factors if f not in extra_primes]

        # Size of 2nd largest prime
        max_prime[len(str(factors[-2:][0]))] += 1
        step += 1

      if not gmpy2.is_prime(t):
        required_steps[1000] += 1
        print ("\t {} Gave({}th time) up on step {}".format(n, required_steps[1000], step))
      else:
        home_primes[(10, n, step-1)] = [t]
        required_steps[step] += 1

  except KeyboardInterrupt:
    print("Stopping from ^C")

  with open("home_primes.txt", "w") as f:
    for base, start, step in sorted(home_primes.keys()):
      factors = home_primes[(base, start, step)]
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
