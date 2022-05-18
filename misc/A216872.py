import gmpy2
from collections import Counter
import tqdm

#s = str(2 ** 880628)
#values = Counter(s).values()
#print ("\t", values, [gmpy2.is_prime(v) for v in values])

HALF = 10

i = 0
#power_two = 2 ** 0
upper = 0
lower = 1

for n in tqdm.tqdm(range(1, 20000)):
    #power_two *= 2
    #t = power_two
    #upper, lower = divmod(t, HALF)

    upper *= 2
    lower *= 2
    if lower > HALF:
        lower -= HALF
        upper += 1

    if upper.bit_length() > lower.bit_length() + 200:
        upper, d = divmod(upper, 10)
        lower += HALF * d
        HALF *= 10

    # 2 ** 877 ~ 10 ** 264
    #if n % 877 == 0:
    #    HALF = 10 ** (264 // 2 * n // 877)
    #    upper, lower = divmod(2 ** n, HALF)

    # Can I make this faster by doing upper half / lower half and handling carry between the two
    # THIS IS SUPER SLOW PART
    values = Counter(str(HALF)).values()
#    values = (Counter(str(upper)) + Counter(str(lower))).values()

    #values = Counter(str(t)).values()
    count_prime = sum(gmpy2.is_prime(v) for v in values)
    if count_prime >= 8:
        values = sorted(values)
        i += 1
        print ("\t", i, n, "\t", values, [v for v in values if not gmpy2.is_prime(v)])
