
with open("b076586_raw.txt") as f:
  data = f.readlines()

with open("b076586.txt") as f:
  bfile = f.readlines()

bfile = dict(list(map(int, l.split())) for l in bfile)
print (bfile)

data = [l.strip() for l in data]

i = 0
while i < len(data):
  line = data[i]
  if not line:
    i += 1
    continue

  timing = ""

  # Timing line
  if "eta" in line:
    timing = line[line.index(" 1:") + 1:]
    i += 1

  # final counts
  assert data[i].startswith("1:"), data[i]
  counts = dict(list(map(int, v.split(":"))) for v in data[i].split())
  i += 1

  if timing:
    for v in timing.split():
      k, c = map(int, v.split(":"))
      assert counts[k] >= c, timing

  n, count = map(int, data[i].split())
  i += 1

  assert sum(counts.values()) == count, (count, counts.values())
  print (n, count)
  assert bfile[n] == count, (n, count, bfile[n])


'''
+  (1928.8m (66007/66148) 99.8%, eta 1933m) 9811523228: 1:25  5:66012  10:14448473  15:243708505  20:764009310  25:682357779  30:223614091  35:31809006  40:2214367  45:82320  50:1701  55:27
+  1:25  2:284  3:2334  4:13971  5:66148  6:258135  7:856474  8:2470221  9:6316224  10:14472680  11:30059972  12:57073262  13:99779261  14:161679487  15:244056172  16:344838699  17:457926639  18:573615153  19:679971388  20:764970905  21:818803633  22:835938468  23:815652908  24:762109592  25:683132309  26:588430064  27:487794742  28:389728773  29:300467548  30:223848129  31:161299926  32:112553410  33:76148712  34:49963760  35:31840122  36:19723008  37:11883408  38:6966585  39:3978450  40:2216366  41:1201557  42:636731  43:330274  44:166761  45:82371  46:39718  47:18942  48:8509  49:3839  50:1701  51:766  52:322  53:137  54:55  55:27  56:8  57:2
+100 9823399067
'''
