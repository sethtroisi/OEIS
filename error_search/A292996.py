import math

temp = math.factorial(220)

def a(n):
    b = {}
    b[1] = 2
    def recursive(k):
        if k in b:
                return b[k]
        r = recursive(k-1)
        z = r * (r + k-1 + n)
        if z % k != 0:
                return 0
        b[k] = (z // k) % temp
        return b[k]

    for k in range(2, 150):
        if recursive(k) == 0:
            return k

s = 0
for n in range(1, 102):
    an = a(n)
    s += an
    print (n, an)
print("\t", s)
