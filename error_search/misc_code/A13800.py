import primesieve
import gmpy2
import time

"""
testing 32 values is 99.999% effective

Can check if there's enough space but 32 value test is faster and more simple.
    Upper bits in S need to miss all lower bits in S
    Calculate a count for both intervals.

    High side [0 ... Am ... p ... max_bit]
         1. Look at how many were set in a previous step
         2. Assume minimum in [p ... max_bit] by having [Am ... p-1] full dense
    Low side [0 ... Am ... max_bit - p ... Am+1]
         1. Look up from previous step
         2. Assume as less dense as possible between [Am ... max_bit - p]

Could possible scan upwards in s looking for a place where two limbs mesh
    See Boyer-Moore
    Would need to skip more than average prime-gap (~10)

    Needs 18 empty from 32, 39 from 64, 54 from 96, 78 from 128

    Can skip max(18 - popcount(limb_a), 39 - popcount(limb_b))
    Probably better to look at byte (or short) and build a jump table

At this point 25% of time is Iterating primes!
"""

# sorted([sum([t*((b>>i)&1) for i, t in enumerate([2,3,7,11,29])]) for b in range(2 ** 5)])
BITS_TO_CHECK = (
    0, 2, 3, 5, 7, 9, 10, 11, 12, 13, 14, 16, 18, 20, 21, 23, 29, 31,
    32, 34, 36, 38, 39, 40, 41, 42, 43, 45, 47, 49, 50, 52, 53, 55, 56, 58, 60, 62, 63
)

BIT_MASK = sum(1 << b for b in BITS_TO_CHECK)
print(f"\tBIT MASK: {BIT_MASK:#b}")

jump_table = []
for i in range(2 ** 16):
    jump = 0
    while i & BIT_MASK:
        i >>= 1
        jump += 1
    jump_table.append(jump)

print(f"\tAvg jump: {sum(jump_table) / 2 ** 16:.2f}")
print(f"\tCount jump=0: {jump_table.count(0)}")
print()

s = gmpy2.mpz(1) # set of sums {0}
set_bits = 1     # number of set bits in s
max_bit  = 0     # largest bit in s

it = primesieve.Iterator()
p = 0
skip_to = 0
for n in range(1, 35+1):
    t = s
    prev_shift = 0
    test_method = [0, 0]
    while True:
        p = it.next_prime()

        if p < skip_to:
            continue

        shift = p - prev_shift

        # Grab 64 bits and do several jumps
        limb = int(t[shift: shift+64])

        # Check if any of of first 39 entries would collide
        if n > 6 and limb & BIT_MASK:
            # limb & BIT_MASK implies a non-zero jump is needed
            jump = max(1, jump_table[limb & 0xFFFF])
            limb >>= jump
            jump2 = jump_table[limb & 0xFFFF]
            limb >>= jump2
            jump3 = jump_table[limb & 0xFFFF]
            skip_to = p + jump + jump2 + jump3

            test_method[0] += 1
            continue

        if p > max_bit:
            break

        # t will likely be very short after this
        t >>= shift
        prev_shift = p

        if test_method[1] == 0:
            print(f"\t first test at {p=}, len(t) = {t.bit_length()}")

        test_method[1] += 1

        # Verify this runs in O(t)
        if (s & t) == 0:
            break

    assert s & (s >> p) == 0

    s += s << p
    set_bits *= 2  # Part of the definition
    max_bit += p

    print(n, p)
    #print(n, p, "\t", max_bit, "\t", test_method)
