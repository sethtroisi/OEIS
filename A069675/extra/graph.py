import matplotlib
import matplotlib.pyplot as plt
import math
import gmpy2
from tqdm import tqdm

results = "results_100000.txt"

with open(results) as results_file:
  numbers = []
  for line in results_file.readlines():
    n = line.split()
    if len(n) == 1:
      numbers.append(int(n[0]))
    else:
      a, b, c = int(n[0]), int(n[4]), int(n[6])
      numbers.append(a * 10 ** b + c)

for i, n in enumerate(tqdm(numbers), 1):
  # Use a big Rabin-Miller confidence
  assert gmpy2.is_prime(n, 100), i


x = list(range(1, len(numbers) + 1))

loglog = [math.log(math.log(n)) for n in numbers]

fig, ax = plt.subplots()
ax.plot(x, loglog, 'k.', markersize=3)
ax.set(xlabel='n', ylabel='ln(ln(a(n)))',
       title='Log Log A069675')
ax.grid()

fig.savefig("LogLogA069675.png")
plt.show()
