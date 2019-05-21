import MathLib
from collections import Counter
from tqdm import tqdm
import multiprocessing

print ("Factoring")
factors = MathLib.sieveOfFactors(10 ** 5)
print ("Factored!")

# Max factors is ~20 (2^20), so number of partition is small (674)

def partition(values):
    if len(values) == 0:
       return

    if len(values) == 1:
        yield ((values[0],),)
        return

    first, *rest = values
    for smaller in partition(rest):
        for n, subset in enumerate(smaller):
            yield smaller[:n] + (((first,) + subset),) + smaller[n+1:]
        yield ((first,),) + smaller

def single(d):
  n, factors = d
  # merge some of the primes together
  for grouping in partition(factors):
#    print ("\tGrouping:", grouping)
    joined = list(map(MathLib.product, grouping))

    counts = Counter(joined)
    modified_sum_of_div = MathLib.product( (p ** (1 + c) - 1) // (p -1 ) for p, c in counts.items() )
    if 2 * n == modified_sum_of_div:
      print (n, "\t", factors, grouping)
      return n
  return None

with multiprocessing.Pool(8) as pool:
  results = pool.map(single, enumerate(factors))
  results = [r for r in results if r is not None]

results = sorted(set(results))
print (results)

