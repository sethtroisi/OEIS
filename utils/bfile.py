import os

def WriteListToFile(number, seq, filename=""):
  # A000001 => 000001
  number = number.lstrip('A')

  filename = filename or "b{}_py.txt".format(number, seq)
  f = open(filename, "w")

  for i, a in enumerate(seq, start=1):
    f.write("{} {}\n".format(i, a))

  f.close()

  size = os.path.getsize(filename)

  human_size = float(size)

  precision = 0
  prefixes = ["bytes", "KB", "MB", "GB"]
  for prefix in prefixes:
    if size < 2 ** 10:
      break
    size /= 2 ** 10
    precision = 1


  print ("\twrote {} lines ({} {}) to {}".format(
      len(seq),
      round(size, precision), prefix,
      filename))
