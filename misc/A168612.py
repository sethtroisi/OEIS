import itertools
import math

an = [1]
for i in range(2, 20):
    for zeros in itertools.count(0):
        binary = "1" + "0" * zeros + bin(an[-1])[2:]
        new_n = int(binary, 2)
        if all(math.gcd(new_n, past) == 1 for past in an):
            an.append(new_n)
            assert binary.endswith(bin(an[-1])[2:])
            print (f"{i} {new_n}")
            break


# Minification attempt

an = [1]
for i in range(2, 20):
    a_i = next(a_i for zeros in itertools.count(0)
            if all(math.gcd(a_i := an[-1] + (1 << (zeros + an[-1].bit_length())),  past) == 1 for past in an))
    an.append(a_i)
    print (f"{i} {a_i}")
