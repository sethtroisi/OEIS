import gmpy2
import math
from multiprocessing import Pool
from tqdm import tqdm

import sys
sys.setrecursionlimit(5000)

def gen_next_gen(test):
    left, right, length = test
    merge = 10 ** length

    new = []
    for x in range(0, 10):
        l_temp = (left * 10 + x)
        join = l_temp * 10 * merge
        for y in range(0, 10):
            r_temp = right + y * merge
            temp = join + r_temp
            if gmpy2.is_prime(temp):
                new.append((l_temp, r_temp, length + 1))
    return new

total_count = 0
global_best = 30
def recursive(test):
    global global_best, total_count

    total_count += 1
    if total_count % 1000000 == 0:
        print("\t", total_count)

    left, right, merge, bitset = test

    if left >= global_best:
        global_best = left
        pairs = bin(bitset).count("1")
        print (2 * pairs, pairs, left, right)

    next_merge = 10 * merge

    new = []
    for x in range(0, 10):
        l_temp = (left * 10 + x)
        join = l_temp * 10 * merge

        for y in range(0, 10):
            new_bitset = bitset | (1 << (10*x+y))
            if bitset == new_bitset:
                continue

            r_temp = right + y * merge
            temp = join + r_temp

            if gmpy2.is_prime(temp):
                r_temp = right + y * merge
                recursive((l_temp, r_temp, next_merge, new_bitset))


# Takes ~200-300 minutes on 10 cores
def puzz_954():
    with Pool() as p:
        gen = [(p//10, p % 10, 1) for p in range(10, 100) if gmpy2.is_prime(p)]
        assert len(gen) == 21, gen

        while gen[0][2] <= 3:
            next_gen = [a for l in p.map(gen_next_gen, gen) for a in l]

            print(gen[0][2], len(next_gen))
            gen = next_gen

#    gen.sort(reverse=True)
    import random
    random.shuffle(gen)

    for test in tqdm(gen):
        # verify pairs are unique here and generate bitset
        left, right, length = test
        strL = str(left)
        strR = str(right).rjust(length, "0")
        assert len(strL) == len(strR)

        pairs = [int(strL[i] + strR[length-i-1]) for i in range(length)]
        if len(pairs) != len(set(pairs)):
            print ("\tskipping", strL + strR, pairs)
            continue

        bitset = 0
        for p in pairs:
            bitset |= 1 << p

        recursive((left, right, 10 ** length, bitset))


puzz_954()
