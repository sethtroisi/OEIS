import gmpy2
import MathLib
import subprocess

def factor_large(n, b1=10**7):
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
  extra_primes = set(map(int, f.readlines()))

# copy to see if we need to write new_factors
original = set(extra_primes)

print ("Found", len(extra_primes), "extra primes")
for p in extra_primes:
  assert gmpy2.is_prime(p), p
print ("Verified")

# Also see A056938
try:
  for n in range(49, 50):
    print (n)
    step = 0
    t = n
    while not gmpy2.is_prime(t):
      s = t
      factors = []
      for p in primes + sorted(extra_primes):
        while t % p == 0:
          factors.append(p)
          t //= p
          if t == 1:
            break

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
      extra_primes.update(factors)
      new = int("".join(map(str, factors)))
      print ("\t", s, factors)
      print ("\t", step, new)
      t = new
except KeyboardInterrupt:
  print("Stopping from ^C")

if extra_primes > original:
  with open("factors", "w") as f:
    for p in sorted(extra_primes):
      f.write(str(p) + "\n")
