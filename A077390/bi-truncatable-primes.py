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
def biTruncatedPrimes(f):
  # Used to store less numbers in interesting.txt
  interesting = 10

  with Pool() as p:
    # Odd length left-and-right-truncated
    iterA = [2,3,5,7]
    assert len(iterA) == 4, iterA

    # Even length left-and-right-truncated
    iterB = [p for p in range(10, 100) if all(p % q != 0 for q in [2,3,5,7])]
    assert len(iterB) == 21, iterB

    total = len(iterA) + len(iterB)

    print_total = 0
    for smallPrimes in [iterA, iterB]:
      for print_total, c in enumerate(smallPrimes, print_total+1):
        f.write("{} {}\n".format(print_total, c))

    assert total == print_total == 25, (total, print_total)

    new_length = 2
    while iterA or iterB:
      new_length += 1

      # Build new primes of new_length digits from iterA

      iterA = [(a, new_length) for a in iterA]

      # Flatten the children lists
      iterC = [a for l in p.map(next_gen, iterA, chunksize=100) for a in l]
      iterC.sort()

      total += len(iterC)

      print(new_length, total, len(iterC), iterC[:3], iterC[-3:])

      for print_total, c in enumerate(iterC, print_total+1):
        if (print_total % interesting == 0 or
            c == iterC[0] or c == iterC[-1] or
            print_total >= 920720000 or print_total <= 200):
          f.write("{} {}\n".format(print_total, c))
          if print_total == interesting * 20:
            interesting = print_total

      iterA, iterB, iterC = iterB, iterC, []


with open("interesting.txt", "w", buffering=1) as f:
  biTruncatedPrimes(f)
