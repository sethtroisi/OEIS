#!/usr/bin/env python3

import urllib.request
import re

URL = "http://chesswanks.com/num/LTPs/"
cache = "html"

#with urllib.request.urlopen(URL) as response:
#  html = response.read()

#with open(cache, "wb") as f:
#  f.write(html)

with open(cache, "rb") as f:
  html = f.read()

html = html.decode()

print (len(html))
#print (html)

#<img src='HTMLFiles/index_89.gif' alt='Graphics:base 89' width='867' height='521' style='vertical-align:middle' /></p>
#<h3 class="Print">{23, 246, 1677, 8357, 33137, 109400, 308223, 756652, 1644822, 3212331, 5689857, 9229972, 13804597, 19151885, 24771356, 30004840, 34189636, 36755721, 37422605, 36169130, 33275833, 29219017, 24530119, 19721482, 15215674, 11286010, 8058216, 5548163, 3686446, 2368648, 1470780, 885767, 517712, 293865, 161266, 86384, 44952, 22865, 11252, 5381, 2494, 1207, 539, 257, 116, 43, 17, 12, 3}</h3>
#<h3 class="Print">largest prime: {24, 70, 76, 6, 26, 2, 18, 78, 38, 72, 54, 72, 24, 14, 48, 78, 42, 78, 30, 46, 60, 86, 26, 86, 12, 54, 20, 62, 86, 20, 14, 44, 30, 24, 42, 78, 66, 28, 30, 2, 62, 42, 22, 58, 84, 84, 78, 84, 83} = 92279250320107991895472932024602610143989437030941920495173615107877916587698486062397943030427</h3>


matches = list(re.finditer(
    r'<img .{100,140}<h3 class="Print">[^<]*</h3><h3 class="Print">[^<]*</h3>',
    html.replace("\n", "")))

print (len(matches))

for match in matches:
  text = match.group()
  base = re.search(":base ([1-9][0-9]*)", text)
  assert base, text
  base = int(base.group(1))

  counts = re.search(">{([0-9, ]+)}", text)
  assert counts, text

  counts = counts.group(1).replace(", ", " ").split()
  assert all(c.isdigit() for c in counts), counts
  counts = list(map(int, counts))
  print (base, sum(counts), len(counts))
