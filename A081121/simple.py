import tqdm

MAX = 10000

cubes = sorted(x*x*x for x in range(0, 10 ** 6 + 2))

found = set()

c_i = 0

for si in tqdm.tqdm(range(0, 10 ** 9)):
  s = si * si

  # start near s
  while cubes[c_i] < s:
    c_i += 1

  for ci in range(c_i, len(cubes)):
    c = cubes[ci]
    n = c - s
    if n  > MAX:
      break
    if n > 0:
      found.add(n)

print ("found {} solutions".format(len(found)))

with open("b081121.txt") as f:
  b_data = [int(line.split()[1]) for line in f.readlines() if line.strip()]

print ("terms: {}, min: {}, max: {}".format(
  len(b_data), min(b_data), max(b_data)))
print ()

not_present = set(range(min(b_data), max(b_data) + 1)) - found
print ("not found:", sorted(set(b_data) - not_present))

extra = not_present - set(b_data)
print ("extra({}): {}".format(len(extra), sorted(extra)))

b_data = sorted(set(b_data) - found)

with open("b081121_new.txt", "w") as f:
  for i, n in enumerate(sorted(b_data), 1):
    f.write("{} {}\n".format(i, n))

