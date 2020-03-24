import gmpy2
import itertools
import subprocess
import math
import time

from collections import defaultdict

from factordb.factordb import FactorDB



START = 2
STOP = 5000


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


def attempt_factorization(s, known_factors):
  t = s
  factors = []

  for factor in known_factors:
    # Last factor maybe non-prime
    if gmpy2.is_prime(factor):
      t //= factor
      factors.append(factor)

  # Toggle to if True: to recheck factordb.
  if t >= 1e10 and t not in known_factors:
    # Check factorDB (probably already been done)
    time.sleep(0.2)
    factordb = FactorDB(t)
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

  min_step = {}
  duplicates = {}

  all_primes = set()
  composites = defaultdict(set)
  for key, factors in home_primes.items():
    for p in factors:
      if gmpy2.is_prime(p):
        all_primes.add(p)
      else:
        composites[key].add(p)

    is_terminal = len(factors) == 1 and factors[0] in all_primes

    s = int("".join(map(str, factors)))
    if s in min_step and not is_terminal:
      # Make sure min step isn't previous step or that's stupid
      if min_step[s] == (key[0], key[1], key[2]-1):
        continue
      duplicates[key] = min_step[s]
    else:
      min_step[s] = key

  print ("Found {} primes, {} composites".format(
      len(all_primes), len(composites)))
  return home_primes, min_step, duplicates, composites


def process(home_primes, composites):
  added = False
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
          print ("Breaking, failed to factor C{}: {}".format(len(str(t)), factordb_format(t)))
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
      else:
        print ("\t {} Gave up on step {}".format(n, step))

  except KeyboardInterrupt:
    print("Stopping from ^C")

  return added


# For use with kernprof -v --line-by-line simple.py
#@profile
def run():
  home_primes, min_step, duplicates, composites = load_from_file()
  added = False

  added = process(home_primes, composites)
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
    print ("|start|step|composite|same as|")
    print ("|-----|----|---------|-------|")
    same = defaultdict(list)
    for key, cfs in composites.items():
      same[tuple(sorted(cfs))].append("HP({}).{}".format(key[1], key[2]))

    merged_count = 0
    for (base, start, step), cfs in composites.items():
      assert (base, start, step+1) not in home_primes
      assert len(cfs) and not gmpy2.is_prime(max(cfs))
      formatted_factors = tuple(factordb_format(c) for c in sorted(cfs))
      key = tuple(sorted(cfs))
      if (base, start, step) not in duplicates:
        same_c = same[key]
        assert same_c[0].startswith("HP({})".format(start)), (key, same_c)
        print ("|HP({})|{}|{}|{}|".format(
            start, step, ", ".join(formatted_factors), " ".join(same_c[1:])))
        merged_count += len(same_c) - 1
      count += 1
    print ("{} numbers ({} merged) <= {} have not yet reached a prime".format(
        count, count - merged_count, STOP))
    print ()
    print ()

  if True:
    print ("### Work")
    print ("---")
    print ()
    # TODO use datetime here
    print ("This is a short list of the smallest (and largest) unfactored numbers as of 2020-03.")
    print ()
    print ("|size|start|step|composite|other factor|")
    print ("|----|-----|----|---------|------------|")
    by_size = sorted((c, key) for key, cfs in composites.items() for c in cfs)
    for c, key in by_size[:30] + by_size[-20:]:
      if key in duplicates:
        continue

      others = home_primes[key][:]
      others.remove(c)
      print ("|c{}|HP({})|step {}|{}|{}|".format(
          len(str(c)), key[1], key[2],
          c,
          " * ".join(map(str, others))))

    print()
    print()

  if True:
    deltas = []

    last = ""
    for (base,start,step),factors in sorted(home_primes.items()):
      assert factors == sorted(factors)
      new = "".join(map(str, factors))
      if step > 1 and (base, start, step) not in duplicates:
        delta = len(new) - len(last)
        deltas.append((delta, int(last), int(new), start, step-1))
      last = new

    # For smallest jump | find biggest number
    # For biggest jumps | find smallest number
    deltas.sort(key=lambda d: (d[0], d[1] if d[0] > 3 else -d[1]))

    print ()
    print ("Home Primes with smallest and largest increase in number of digits")
    print ()
    print ("|+digits|HP|current|next|link|")
    print ("|-------|--|-------|----|----|")
    for delta, s1, s2, start, step in deltas[:15] + deltas[-15:]:
      print("|{}|{}|{}|{}|{}|".format(
        delta,
        f"HP({start}).{step}",
        factordb_format(abs(s1)),
        factordb_format(abs(s2)),
        "[FactorDB](http://factordb.com/aliquot.php?type=10&aq={}&big=1)".format(start)))

run()
