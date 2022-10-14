import array
import math
import time

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


def iterate(START, STOP, SEGMENT):
    for start in range(START, STOP, SEGMENT):
        sigmas = SegmentedSieveOfSigma(start, SEGMENT+2)
        for i in range(0, SEGMENT):

            if sigmas[i+2] - sigmas[i] == 2:
                print_match(start + i + 1)

                # Weird non prime terms
                mid = i + 1
                if mid in (435, 8576, 8826):
                    continue

                assert sigmas[i+2] == (start+i+3), (mid+1, sigmas[i+2])
                assert sigmas[i] == (start+i+1),   (mid-1, sigmas[i])

START = 0
SEGMENT = 2 ** 17
STOP = 10 ** 9
iterate(START, STOP, SEGMENT)
