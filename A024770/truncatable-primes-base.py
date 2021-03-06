import gmpy2
from multiprocessing import Pool
from functools import partial

# Left- or Right- truncatable primes
LEFT = True
LEFT = False


SMALL_PRIMES = [2]
while SMALL_PRIMES[-1] < 1000:
  SMALL_PRIMES.append(int(gmpy2.next_prime(SMALL_PRIMES[-1])))


def next_right(cur, base):
  # Add digits to the right
  test = cur * base

  new = []
  for right in range(1, base, 2 - (base % 2)):
    temp = test + right
    if gmpy2.is_prime(temp):
      new.append(temp)
  return new


def next_left(cur, base, new_length):
  # Add digits to the left
  mul = base ** (new_length - 1)

  new = []
  for left in range(1, base):
    temp = left * mul + cur
    if gmpy2.is_prime(temp):
      new.append(temp)
  return new


def truncated_primes(base):
  with Pool() as pool:
    iter_a = [p for p in SMALL_PRIMES if p < base]
    assert SMALL_PRIMES[-1] > base

    total = len(iter_a)

    new_length = 1
    while iter_a:
      new_length += 1

      if LEFT:
        next_left_fun = partial(next_left, base=base, new_length=new_length)
        next_gen = pool.map(next_left_fun, iter_a, chunksize=100)
      else:
        next_right_fun = partial(next_right, base=base)
        next_gen = pool.map(next_right_fun, iter_a, chunksize=100)

      # Flatten the children lists
      iter_a = [a for l in next_gen for a in l]
      iter_a.sort()

      total += len(iter_a)

#      print(new_length, total, len(iter_a), iter_a[:3], iter_a[-3:])

  print (base, total)


########


for base in range(2, 60 + 1):
  # Don't save any of the results
  truncated_primes(base, None)
