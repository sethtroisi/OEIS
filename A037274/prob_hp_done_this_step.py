import math
import gmpy2
from tqdm import tqdm

SIMULATIONS = 5000

#prefix = "71741"
# c152 from HP(3466)
#cXXX = 25695243121353635259858419692719847767472012810564789500454217744414329594815974108332158170426874886245124487058689023388754418742928213371835049836031

prefix = "7316249"
# c251 from HP(49)
cXXX = 26633090926792263436736904630531520479768742849435097127754634822168395438250791486509180275478812779959346908131589660697709489852830934711978704681639399323263270697821325581691729538877317736626598036703631970679737664887720652086830617767029002763

small_prime = 30

assert cXXX > 10 ** (2 * small_prime + 10)
pXX = gmpy2.next_prime(10 ** small_prime)

is_prime = 0
is_prime_theory = 0
for _ in tqdm(range(SIMULATIONS)):
  # Simulate factoring cXXX into pXX * pOther

  pXX = int(1.001 * pXX)
  pOther = cXXX // pXX

  # Make primes have prime digits
  pXX    = gmpy2.next_prime(pXX)
  pOther = gmpy2.next_prime(pOther)

  step = int(prefix + str(pXX) + str(pOther))
#  assert 258 <= len(step) <= 261, (len(step), len(pXX))
  if gmpy2.is_prime(step):
    is_prime += 1
  is_prime_theory += 1 / math.log(step)


print ("Finish after this step in {}/{} = {:.2%} of the simulations".format(
    is_prime, SIMULATIONS, is_prime / SIMULATIONS))
print ("Using is_prime = 1 / log(n): {:.2%}".format(
    is_prime_theory / SIMULATIONS))
