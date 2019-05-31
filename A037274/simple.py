import gmpy2
import itertools
import subprocess
import math
import time


from factordb.factordb import FactorDB

from collections import Counter, defaultdict

# Also see A056938

def product(factors):
  temp = 1
  for factor in factors:
    temp *= factor
  return temp


def factordb_format(number):
  if number < 1e10:
    return str(number)
  strN = str(number)
  length = len(strN)
  if number < 1e24:
    return "{}<{}>".format(strN, length)
  return "{}...{}<{}>".format(strN[:10], strN[-2:], length)

def split_to_lines(number, max_size):
  size = max_size - 2

  # split this number evenly over multiple lines
  needed_lines = (len(number) - 1) // size + 1
  assert size * needed_lines >= len(number)

  # split evenly onto this many lines
  per_line = len(number) // needed_lines
  # this many lines get 1 extra
  extra = len(number) % needed_lines
  assert per_line + (extra > 0) <= size

  lines = []
  for l in range(1, needed_lines+1):
    # take per_line, plus potentially one extra
    this_line = number[:per_line + (extra > 0)]
    number = number[len(this_line):]
    this_line += " /" if l != needed_lines else ""
    lines.append(this_line)

  return lines

def row_format(string, max_size=60):
  if len(string) <= max_size:
    return string

  mult = " * "
  if mult in string:
    parts = string.split(mult)

    lines = []
    line = ""
    for part in parts:
      merged = line + part + mult
      if len(merged) <= max_size + 1: # trailing space
        line = merged
        continue
      elif line:
        lines.append(line.strip())
        line = ""

      assert line == ""
      if len(part) <= max_size - 2:
        lines.append(part + " *")
        continue

      lines += split_to_lines(part + " *", max_size)

    temp = "<br>".join(lines)
    assert temp.endswith(" *"), temp[-20:]
    return temp[:-2]

  return "<br>".join(split_to_lines(string, max_size))


def factor_large(n, b1=10**6):
  args = ["ecm", "-q", "-c", "10", str(b1)]
  print ("\t\t", " ".join(args))
  result = subprocess.run(
      args,
      input=str(n).encode(),
      stdout=subprocess.PIPE)

  if result.returncode == 8:
    # Need to rerun with smaller b1
    print("\t\tfound self ({} with b1={})".format(n, b1))
    return factor_large(n, b1= max(100, b1 // 90))

  return list(map(int, result.stdout.strip().split()))


def attempt_factorization(s, know_factors):
  t = s
  factors = []

  for factor in know_factors:
    # Last factor maybe non-prime
    if gmpy2.is_prime(factor):
      t //= factor
      factors.append(factor)

  #if False:
  if t >= 1e10 and t not in know_factors:
    # Check factorDB (probably already been done)
    time.sleep(0.2)
    factordb = FactorDB(s)
    factordb.connect()
    factordb_factors = factordb.get_factor_list()
    if factordb_factors:
      print ("\t\tfactordb:", factordb.get_status(), factordb_factors)
      for factor in factordb_factors:
        if gmpy2.is_prime(factor):
          t //= factor
          factors.append(factor)

  # small trial division
  p = 2
  while t > 1 and t < 1e10:
    while t % p == 0:
      t //= p
      factors.append(p)
      if t == 1:
        break
    p += 1 + (p&1)

  return t, factors


def load_from_file():
  home_primes = defaultdict(list)

  n = None
  s = None
  with open("home_primes.txt") as f:
    # each line is "<base> <start> <step> <status>: <factor> <factor> ..."
    for line in f.readlines():
      pre, post = line.strip().split(":")
      *pre, status = pre.split()
      base, start, step, = map(int, pre)

      if start != n:
        n = start
        s = n

      factors = list(map(int, post.split()))
      assert status in ("FF", "P", "CF"), line
      home_primes[(base, start, step)] = factors
      assert product(factors) == s, (start, step, s, factors)
      s = int("".join(map(str, factors)))

  all_primes = set()
  composites = defaultdict(set)
  for key, factors in home_primes.items():
    for p in factors:
      if gmpy2.is_prime(p):
        all_primes.add(p)
      else:
        composites[key].add(p)

  print ("Found {} primes, {} composites".format(
      len(all_primes), len(composites)))
  return home_primes, all_primes, composites


# For use with kernprof -v --line-by-line simple.py
#@profile
def run():
  START = 2
  STOP = 5000

  home_primes, all_primes, composites = load_from_file()
  added = False

  required_steps = Counter()
  try:
    for n in range(START, STOP+1):
      print (n)
      t = n
      for step in itertools.count(1):
        if gmpy2.is_prime(t):
          break

        s = t

        key = (10, n, step)
        original = home_primes[key]
        t, factors = attempt_factorization(s, original)

        factors.sort()

        if t > 1:
          # t is composite
          factors.append(t)
          composites[key].add(t)

        assert product(factors) == s, (s, t, factors)

        if factors != original:
          home_primes[key] = factors
          added = True
          print ("\t\tnew factor", factors)

        if t > 1:
          print ("Breaking, failed to factor C{}: {}".format(len(str(t)), t))
          break

        new = int("".join(map(str, factors)))
        t = new

        if False:
            if gmpy2.is_prime(s):
              if new < 1e40:
                print ("\t", step, new, "from", s, factors)
              else:
                print ("\t", step, new)
                print ("\t\tfrom", factors)

      if gmpy2.is_prime(t):
        home_primes[(10, n, step)] = [t]
        required_steps[step] += 1
      else:
        required_steps[1000] += 1
        print ("\t {} Gave({}th time) up on step {}".format(n, required_steps[1000], step))

  except KeyboardInterrupt:
    print("Stopping from ^C")

  if added:
    with open("home_primes.txt", "w") as f:
      for base, start, step in sorted(home_primes.keys()):
        factors = home_primes[(base, start, step)]
        if not factors:
          continue

        if all(gmpy2.is_prime(f) for f in factors):
          if len(factors) == 1:
            status = "P"
          else:
            status = "FF"
        else:
          status = "CF"

        f.write("{} {} {} {}: {}\n".format(
            base, start, step, status, " ".join(map(str, factors))))

  # Sections copied into README.md
  if True:
    ranges = [(2,100), (2,499)] + [(a*500, a*500 + 499) for a in range(1, STOP//500)]
    for low, high in ranges:
      filename = "RESULTS_{}_{}.md".format(low, high)
      print ("Genarating", filename)
      template = """
## [Back](../README.md)

## Results for A037274 a(n) n={}..{}
---

|start|step|number|factors|
|-----|----|------|-------|
{}

"""
      rows = []
      for (_,start,step),factors in sorted(home_primes.items()):
        if start not in range(low, high+1):
          continue

        num = row_format(str(product(factors)), max_size=40)
        if len(factors) == 1:
          factors = "Home Prime!" if gmpy2.is_prime(min(factors)) else "Unfactored composite"
        else:
          mult = " * ".join(map(str, sorted(factors)))
          factors = row_format(mult, max_size=50)

        columns = [start, step, num, factors]
        rows.append("|" + "|".join(map(str, columns)) + "|")

      with open("results/" + filename, "w") as f:
        f.write(template.format(
            low, high,
            "\n".join(rows)))

  if True:
    count = 0
    print ()
    print ()
    print ("### Unterminated")
    print ("---")
    print ()
    # Move the "These <X> a(n) that have not..." line here
    print ()
    print ("|start|step|composite|")
    print ("|----------|----|---------|")
    for (base, start, step), cs in sorted(home_primes.items()):
      if (base, start, step+1) in home_primes:
        continue
      if not (len(cs) == 1 and gmpy2.is_prime(max(cs))):
        cfs = [factordb_format(c) for c in cs if not gmpy2.is_prime(c)]
        print("|HP({})|{}|{}|".format(
          start, step, ", ".join(cfs)))
        count += 1
    print ("These", count, "a(n) (n = {}..{}) have not yet reached a prime".format(
        START, STOP))
    print()
    print()

  if True:
    print ("### Work")
    print ("---")
    print ()
    print ("This is a short list of the smallest (and largest) unfactored numbers as of 2019-05.")
    print ()
    print ("|size|start|step|composite|other factor|")
    print ("|----|-----|----|---------|------------|")
    by_size = sorted((c, key) for key, cs in composites.items() for c in cs)
    for c, key in by_size[:30] + by_size[-20:]:
      others = home_primes[key]
      others.remove(c)
      print ("|c{}|HP({})|step {}|{}|{}|".format(
          len(str(c)), key[1], key[2],
          c,
          " * ".join(map(str, others))))
    print()
    print()

run()
