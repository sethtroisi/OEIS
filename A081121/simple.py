MAX = 10000

squares = sorted(y*y for y in range(0, 33000))
cubes = sorted(x*x*x for x in range(-1000, 1000))

found = set()

for s in squares:
  # start near s
  for c in cubes:
    n = c - s
    if n  > MAX:
      break
    if n > 0:
      found.add(n)

print ("found {} solutions".format(len(found)))

with open("b081121.txt") as f:
  b_data = [int(line.split()[1]) for line in f.readlines() if line.strip()]

print (len(b_data), min(b_data), max(b_data))

to_remove = []
for n in b_data:
  if n in found:
    print ("What what", n)
    to_remove.append(n)

b_data = [n for n in b_data if n not in to_remove]

with open("b081121_new.txt", "w") as f:
  for i, n in enumerate(sorted(b_data), 1):
    f.write("{} {}\n".format(i, n))

