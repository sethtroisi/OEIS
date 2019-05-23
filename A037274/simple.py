import gmpy2
import MathLib
import subprocess

from collections import Counter

# Also see A056938

def factor_large(n, b1=10**6):
  args = ["ecm", "-q", "-I", "1000", "-c", "10", str(b1)]
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

primes = MathLib.sieveOfErat(10 ** 5)
with open("factors") as f:
  extra_primes = sorted(map(int, f.readlines()))

# copy to see if we need to write new_factors
original = set(extra_primes)

print ("Found", len(extra_primes), "extra primes")
for p in extra_primes:
  assert gmpy2.is_prime(p), p
print ("Verified")

required_steps = Counter()
max_prime = Counter()
try:
  for n in range(2, 10000):
    print (n)
    step = 0
    t = n
    while not gmpy2.is_prime(t) and len(str(t)) < 30:
      s = t
      factors = []
      for p in extra_primes:
        while t % p == 0:
          factors.append(p)
          t //= p
          if t == 1 or p > t:
            break

      for p in primes:
        while t % p == 0:
          factors.append(p)
          t //= p
          if t == 1 or p*p > t:
            break

      if t > 1:
        print ("\t\tpre:", factors, t)
      while t > 1:
        if gmpy2.is_prime(t):
          factors.append(t)
          t //= t
          continue

        ecm_factors = sorted(factor_large(t))
        print ("\t\tecm:", ecm_factors)
        for f in ecm_factors:
          if gmpy2.is_prime(f):
            t //= f
            factors.append(f)

      step += 1
      factors.sort()
      new = int("".join(map(str, factors)))
      print ("\t", step, new, "from", s, factors)
      t = new

      extra_primes += [f for f in factors if f not in extra_primes]
      extra_primes.sort()

      # Size of 2nd largest prime
      max_prime[len(str(factors[-2:][0]))] += 1

    if not gmpy2.is_prime(t):
      print ("\t Gave up on step", step)
      required_steps[0] += 1
    else:
      required_steps[step] += 1

except KeyboardInterrupt:
  print("Stopping from ^C")

if set(extra_primes) > original:
  with open("factors", "w") as f:
    for p in sorted(extra_primes):
      f.write(str(p) + "\n")

print ()
print ("0 steps = stopped")
for s, c in sorted(required_steps.items()):
  print ("{} with {} steps".format(c, s))
print ()
for p, c in sorted(max_prime.items()):
  print ("2nd largest primes had {} digits x {} steps".format(p, c))

