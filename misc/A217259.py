import array
import itertools
import math
import time
from multiprocessing import Pool

from sympy.ntheory.primetest import isprime
from sympy.ntheory import divisor_sigma

def SegmentedSieveOfSigma(start, N):
    """Calculate sigma[i] for i in range(start, start+n)."""

    past = start + N
    sigmas = array.array('l', [1]) * N

    if start == 0:
        sigmas[0] = 0
        sigmas[1] = 0

    for factor in range(2, math.isqrt(past-1) + 1):
        f2 = factor * factor
        assert f2 < past

        if f2 >= start:
            assert f2 - start < N
            # n=factor^2 only gets to count factor, not factor + n/factor
            sigmas[f2 - start] += factor
            ceil = factor + 1
            next_index = f2 + factor - start
        else:
            ceil = (start-1) // factor + 1
            next_index = ceil * factor - start

        for count, index in enumerate(range(next_index, N, factor), ceil):
            # count = number / factor
            sigmas[index] += factor + count

    # Adjust to include n as a divisor of n
    for i in range(N):
        sigmas[i] += start + i

    return sigmas


S = time.time()
found = 0
print_mult = 1
next_time = S + 2

def print_match(mid):
    global found, print_mult, next_time
    found += 1
    if found % print_mult == 0:
        print(f"{found:<8} {mid:<13,}")
        if found == 10 * print_mult:
            print_mult *= 10

    elif time.time() > next_time:
        elapsed = time.time() - S
        rate = mid / elapsed / 10 ** 6
        print(f"{found:<8} {mid:<13,}\t{elapsed:.1f} elapsed {rate:.1f}M/s")
        next_time += 5

    if mid in (435, 8576, 8826):
        return

    assert isprime(mid - 1), mid
    assert isprime(mid + 1), mid

def iterate(START, STOP, SEGMENT):
    # sigma "-2, -1"
    if START <=  1:
        last_sigmas = [-100, -100]
    elif START == 1:
        last_sigmas = [-100, divisor_sigma(1)]
    else:
        last_sigmas = [divisor_sigma(START-2), divisor_sigma(START-1)]

    for start in itertools.count(START, SEGMENT):
        sigmas = SegmentedSieveOfSigma(start, SEGMENT)

        if start > STOP:
            break

        if sigmas[0] - last_sigmas[0] == 2:
            print_match(start - 1)

        if sigmas[1] - last_sigmas[1] == 2:
            print_match(start)

        for i in range(1, SEGMENT-1):
            if sigmas[i+1] - sigmas[i-1] == 2:
                print_match(start + i)

        last_sigmas = list(sigmas[-2:])
        del sigmas

START = 0
SEGMENT = 2 ** 17
STOP = 10 ** 15
iterate(START, STOP, SEGMENT)

# N < 10^8 | 5.3 seconds w/ SEGMENT = 2^16 | 440315 99999588
# N < 10^9 | 48  seconds w/ SEGMENT = 2^16 | 3424509 999999192
# N < 10^9 | 56  seconds w/ SEGMENT = 2^20 | 3424509 999999192
