import primesieve
import gmpy2
import time

s = gmpy2.mpz(1) # set of sums {0}
set_bits = 1     # number of set bits in s
max_bit  = 0     # largest bit in s

it = primesieve.Iterator()
p = 0
for n in range(1, 35+1):
    t = s
    prev = 0
    while True:
        p = it.next_prime()

        # Check the first 16 entries (lets us skip very often)
        # sorted([sum([t*((b>>i)&1) for i, t in enumerate([2,3,7,11])]) for b in range(2 ** 4)])
        # TODO: very speedy version would use a bit mask here
        if n > 4 and any(t.bit_test(i + (p-prev)) for i in (0,2,3,5,7,9,10,11,12,13,14,16,18,21,23)):
            continue

        # TODO could shorten s each time too
        t >>= p - prev
        prev = p

        # Would be nice for this to run faster than O(s)
        if (s & t) == 0:
            break

    '''
    while True:
        p = gmpy2.next_prime(p)

        if p > max_bit:
            break

        # Check if even theoretically possible
        space = max_bit - p

        # Upper bits in S need to miss all lower bits in S
        # Calculate a count for both intervals.
        # High side [0 ... Am ... p ... max_bit]
        #     1. Look at how many were set in a previous step
        #     2. Assume minimum in [p ... max_bit] by having [Am ... p-1] full dense
        # Low side [0 ... Am ... max_bit - p ... Am+1]
        #     1. Look up from previous step
        #     2. Assume as less dense as possible between [Am ... max_bit - p]

        # High side
        #   1. | Could cache this
        previous_a, set_previous_a = max((mb, sb) for mb, sb in An if mb <= p)
        #   2. | All values above p are dense (could recursively
        #   Max possible number of entries < p  =>  Min possible entry >= p
        remaining_top = max(1, set_bits - ((p - previous_a - 1) + set_previous_a))

        previous_b, set_previous_b = max((mb, sb) for mb, sb in An if mb <= space)
        # TODO can add a few more values by using next value and doing density
        remaining_bottom = set_previous_b

        combined = remaining_top + remaining_bottom

        # Now check if possible to have no overlap
        print("\t", p, f"\t {set_previous_a} <= {previous_a}  =>  {remaining_top} >= {p} | {set_previous_b} <= {previous_b} | {combined} vs {space}")

        if combined >= space:
            continue

        if s & (s >> p) == 0:
            break
    '''

    assert s & (s >> p) == 0

    s += s << p
    set_bits *= 2  # Part of the definition
    max_bit += p

    print(n, p, "\t", max_bit)
