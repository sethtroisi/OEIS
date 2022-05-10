import sys
sys.path.append("../utils")
import bfile

STOP = 20000

seen = set([0])

for a in range(1, STOP):
  c = a ** 3

  if c > STOP:
    break

  to_update = []
  for s in seen:
    t = s + c
    if t > STOP:
      continue

    if t not in seen:
      to_update.append(t)

  seen.symmetric_difference_update(to_update)
  print (f"{a=} {len(seen)=}")

seq = sorted(set(range(1, STOP + 1)) - seen)
#for i, s in enumerate(seq, 1):
#  print (i, s)

print (len(seen), len(seq))

bfile.WriteListToFile('A001476', seq, "b001476.txt") 
