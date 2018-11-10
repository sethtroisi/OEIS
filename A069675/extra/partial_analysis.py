import re
from collections import defaultdict

filter_file = "filter_1_100000.partial"

abToD = defaultdict(lambda: [[], []])

for line in open(filter_file).readlines():
  match = re.match('^([0-9]*), ([0-9]*), ([0-9]*): ([0-9-]*)$', line)
  if match:
    d, a, b, res = map(int, match.groups())

    if d <= 5:
      continue

    if res == -1:
      abToD[(a,b)][0].append(d)
    elif res == -2:
      abToD[(a,b)][1].append(d)
  else:
    assert False, line


def filterfalse(a, b):
  return filter(lambda c : not(a(c)), b)


# Function should return FALSE when it can't be prime.
def testHypoth(function, true_set, false_set):
  if len(true_set) == 0:
    assert False

  confirm = len(list(filter(function, true_set)))
  deny    = len(list(filterfalse(function, true_set)))
  usefull = len(list(filter(function, false_set)))
  useless = len(list(filterfalse(function, false_set)))

  if confirm == 0:
    # swap hypothesis around
    confirm, deny = deny, confirm
    usefull, useless = useless, usefull

  # Can't overfilter the true set.
  return deny == 0, (usefull, useless)



for (a, b), (prime, composite) in abToD.items():
#  print (a, b)
#  if (a, b) in ((4, 7), (9, 7)):
#    print ("\t", testHypoth(
#        lambda d : (d % 2 == 0) and (((a * pow(10, d, 28) + b) % 28) not in (1,9,11,15,23,25)),
#        prime,
#        composite))
#
  if (a, b) in ((5, 9),):
    # 5 * x^2 + y^2
    print ("\t", testHypoth(
        lambda d : (d % 2 == 0) and (((a * pow(10, d, 20) + b) % 20) not in (1,9)),
        prime,
        composite))

  if (a, b) in ((8, 7),):
    # 2x^2 + 7y^2
    print ("\t", testHypoth(
        lambda d : (d % 2 == 0) and (((a * pow(10, d, 20) + b) % 20) not in (1,9)),
        prime,
        composite))
