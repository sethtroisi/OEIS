from fractions import Fraction

START_P = 100
POWER_CUTOFF = 100

def frac_m(p, e):
  pe = p ** e
  return Fraction(pe * p - 1, (p-1) * pe)

# Smallest number can go first?
candidates = [(Fraction(1, 1), 1, tuple())]
generation = 0
while generation <= 6:
  generation += 1
  print(generation, len(candidates))
  next_gen = []
  for frac, p, n in candidates:

#    num = frac.numerator
#    for m in range(num, 100*num+1, num):
    for m in range(p+1, START_P):
      #if m < p: continue
      if m % 2 == 0: continue

      for e in range(1, 2 + 1 * (m <= POWER_CUTOFF)):
        me = m ** e
        new_frac = frac * frac_m(m, e)
        new_n = n + ((m,e),)

        if new_frac == 2:
          print("\t",new_n)
          continue

        if new_frac.numerator + 1 == 2 * new_frac.denominator:
          print("\t", new_frac, new_n)

        if new_frac < 2:
          next_gen.append((new_frac, m, new_n))

  candidates = next_gen

'''

5/4 * 6/5 * 7/6 * 8/7 = 30/20 * 8/6 = 3/2 * 8/6 = 24/12 = 2

4 => 5/4

[5?] 10? 15?

5 => 5/4 * 6/5 = 6/4 = 3/2

3? [6?] 9?

6 => 3/2 * 7/6 = 7/4

[7?] 14? 21?

7 => 7/4 * 8/7 = 8/4 = 2/1 !!!
4*5*6*7 = 840

----

descart

3^2  => ((3^3 - 1)/2)/3   = 13/3
7^2  => ((7^3 - 1)/6)/7   = 57/7
11^2 => ((11^3 -1)/10)/11 = 133/11
13^2 => ((13^3 -1)/12)/13 = 186/13


'''
