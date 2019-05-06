import gmpy2
from multiprocessing import Pool

def right_gen(c):
  new = []
  temp = 10 * gmpy2.mpz(c)
  for right in [1,3,7,9]:
    test = temp + right
    if gmpy2.is_prime(test):
      new.append(test)
  return new

def next_gen(test):
  c, mul = test
  new = []
  for left in range(1, 10):
    l_temp = left * mul + 10 * c
    for right in [1,3,7,9]:
      temp = l_temp + right
      if gmpy2.is_prime(temp):
        new.append(temp)
  return new


def rightTruncatedPrimes():
  with Pool(10) as p:
    current = [2,3,5,7]
    for iteration in range(2, 10):
      newGen = p.map(right_gen, current, 10)
      current = [a for n in newGen for a in n]

      print(iteration, len(current), min(current, default=0))

#rightTruncatedPrimes()

def biTruncatedPrimes(f):
  with Pool() as p:
    current = [2,3,5,7]
    assert len(current) == 4, current

    current = [p for p in range(10, 100) if all(p % q != 0 for q in [2,3,5,7])]
    assert len(current) == 21, current

    total = len(current)

    iteration = len(str(current[0]))
    while len(current) > 0:
      iteration += 2

      mul = 10 ** (iteration - 1)
      current = [(c, mul) for c in current]
      newGen = p.map(next_gen, current, 10)
      current = [a for n in newGen for a in n]
      current.sort()

      for c in current:
        f.write(str(c) + "\n")

      total += len(current)
      print(iteration, total, len(current), current[:5], current[-5:])

  # Takes ~100-300 minutes on 10 cores

with open("bi-1.txt", "w") as f:
  biTruncatedPrimes(f)

