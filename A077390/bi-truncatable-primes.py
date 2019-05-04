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
  with Pool(10) as p:
    #current = [2,3,5,7]
    current = [p for p in range(10, 100) if all(p % q != 0 for q in [2,3,5,7])]
    assert len(current) == 21, current

    total = len(current)

    iteration = len(str(current[0]))
    for iteration in range(iteration + 2, 120, 2):
      mul = 10 ** (iteration - 1)
      current = [(c, mul) for c in current]
      newGen = p.map(next_gen, current, 10)
      current = [a for n in newGen for a in n]
      current.sort()

      for c in current:
        f.write(str(c) + "\n")

      total += len(current)
      print()
      print(iteration, total, len(current), current[:5], current[-5:])

  # Takes ~100-300 minutes on 10 cores

with open("bi-2.txt", "w") as f:
  biTruncatedPrimes(f)


  """
  3 59
  5 494
  7 2833
  9 12718
  11 46807
  13 141900
  15 377342
  17 877916
  19 1783865
  21 3332745
  23 5693162
  25 8840023
  27 12827559
  29 17284717
  31 21660200
  33 25647604
  35 28554805
  37 29673844
  39 29660754
  41 28530978
  43 25687627
  45 22285448
  47 18647299
  49 14610917
  51 11147972
  53 8208662
  55 5697725
  57 3911948
  59 2609130
  61 1639381
  63 1015824
  65 604159
  67 345710
  69 196217
  71 107676
  73 56658
  75 29829
  77 15099
  79 7129
  81 3404
  83 1503
  85 688
  87 312
  89 132
  91 46
  93 24
  95 13
  97 3
  99 0
  """


