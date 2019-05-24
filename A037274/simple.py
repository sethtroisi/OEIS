import gmpy2
import subprocess
import math

from collections import Counter

# Also see A056938

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


# For use with kernprof -v --line-by-line simple.py
#@profile
def run():
  with open("factors") as f:
    extra_primes = sorted(map(int, f.readlines()))

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
    for n in range(2, 10000):
      print (n)
      step = 0
      t = n
      while not gmpy2.is_prime(t):
        s = t

        # This is the slowest part of rerunning
        t, factors = trial_factor(t, primes, original)

        if len(str(t)) > 60: break

        if t > 1:
          print ("\t\tpre:", factors, t)

          while t > 1:
            if gmpy2.is_prime(t):
              factors.append(int(t))
              t //= t
              break

            ecm_factors = sorted(factor_large(t))
            print ("\t\tecm:", ecm_factors)
            for f in ecm_factors:
              if gmpy2.is_prime(f):
                t //= f
                factors.append(int(f))


        step += 1
        factors.sort()
        new = int("".join(map(str, factors)))
        print ("\t", step, new, "from", s, factors)
        # For reporting to factor db.
        # for f in factors:
        #   if len(str(s)) > 60 and len(str(f)) > 20:
        #     print (s, "=", f)

        t = new

        extra_primes += [f for f in factors if f not in extra_primes]

        # Size of 2nd largest prime
        max_prime[len(str(factors[-2:][0]))] += 1

      if not gmpy2.is_prime(t):
        required_steps[1000] += 1
        print ("\t {} Gave({}) up on step {}".format(n, required_steps[0], step))
      else:
        required_steps[step] += 1

  except KeyboardInterrupt:
    print("Stopping from ^C")

  if set(extra_primes) > original:
    print("Adding", len(extra_primes) - len(original), "new primes")
    with open("factors", "w") as f:
      for p in sorted(extra_primes):
        f.write(str(p) + "\n")

  print ()
  for s, c in sorted(required_steps.most_common()):
    print ("{} with {} steps".format(c, s))
  print ("1000 steps = stopped early")

  print ()
  for p, c in sorted(max_prime.items()):
    print ("2nd largest primes had {} digits in {} steps".format(p, c))


run()
