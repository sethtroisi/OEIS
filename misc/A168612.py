import itertools
import math

an = [1]
for i in range(2, 100):
    for zeros in itertools.count(0):
        binary = "1" + "0" * zeros + bin(an[-1])[2:]
        new_n = int(binary, 2)
        if all(math.gcd(new_n, past) == 1 for past in an):
            an.append(new_n)
            assert binary.endswith(bin(an[-1])[2:])
            print (f"{i} {new_n}")
            break
