import gmpy2
from multiprocessing import Pool

def next_gen(test):
  c, new_length = test
  mul = 10 ** (new_length - 1)

  new = []
  for left in range(1, 10):
    l_temp = left * mul + 10 * c
    for right in [1,3,7,9]:
      temp = l_temp + right
      if gmpy2.is_prime(temp):
        new.append(temp)
  return new


# Takes ~200-300 minutes on 10 cores
def biTruncatedPrimes(f, f2):
  # Used to store less numbers in interesting.txt
  interesting = 10

  with Pool() as p:
    # Odd length left-and-right-truncated
    iterA = [2,3,5,7]
    assert len(iterA) == 4, iterA

    # Even length left-and-right-truncated
    iterB = [p for p in range(10, 100) if all(p % q != 0 for q in [2,3,5,7])]
    assert len(iterB) == 21, iterB

    total = 1
    for total, c in enumerate(iterA, total):
      f.write(str(c) + "\n")
      f2.write("{} {}\n".format(total, c))
    for total, c in enumerate(iterB, total+1):
      f.write(str(c) + "\n")
      f2.write("{} {}\n".format(total, c))

    assert total == 25, total

    new_length = len(str(iterA[0])) + 1
    while iterA or iterB:
      new_length += 1

      iterA = [(a, new_length) for a in iterA]

      # Flatten
      iterC = [a for l in p.map(next_gen, iterA, chunksize=100) for a in l]
      iterC.sort()

      print(new_length, total + len(iterC), len(iterC), iterC[:3], iterC[-3:])

      for total, c in enumerate(iterC, total+1):
        if total % interesting == 0 or c == iterC[0] or c == iterC[-1]:
          f2.write("{} {}\n".format(total, c))
          if total == interesting * 10:
            interesting = total

        f.write(str(c) + "\n")

      iterA, iterB, iterC = iterB, iterC, []


with open("left-right-truncated.txt", "w") as f, open("interesting.txt", "w") as f2:
  biTruncatedPrimes(f, f2)

