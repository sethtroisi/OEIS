import array
import math
import sys
import time


def get_prime_array(n):
    """Get prime array

    use primesieve if python, use sympy if pypy
    """
    try:
        import primesieve
        return primesieve.primes(n)
    except:
        from sympy import sieve
        a = array.array('Q')
        a.extend(sieve.primerange(2, n+1))
        return a


def A000047(bits: int) -> int:
    # range = [0, last), manually correct for 2 ** n at the end
    last  = 2 ** bits

    primes = get_prime_array(last)
    #print(f"Primes {primes[0]} to {primes[-1]}")

    # easier than array
    status = bytearray(last)
    for i in range(last):
        status[i] = 1

    max_e = n
    for p in primes:
        if p % 8 not in (3, 5):
            continue

        p2 = p * p
        pp = p
        for e in range(1, max_e+1, 2):
            if pp > last:
                max_e = e
                break

            assert pp % 8 in (3, 5)
            # Mark off all multiples (not of p)
            # Only need to mark off multiples that are valid

            # Only ~1/3 are valid but takes longer doing multiplication
            #for i in range(1, last // pp + 1):
            #    if status[i] == 1 and i % p != 0:
            #        status[i * pp] = 0

            status[pp] = 0
            m = 0
            while m < last:
                # Mark off (i * p + j) * pp, 1 <= j < p

                # Skip multiple p * pp
                stop = min(last, m + pp * p)
                m += pp
                for m in range(m, stop, pp):
                    status[m] = 0
                m += pp

            pp *= p2

    # -1 for 0, +1 for 2^n
    count = sum(status)
    #for i, s in enumerate(status):
    #    if s: print(i)
    return count


def A000047_fast(bits: int) -> int:
    n = 2 ** bits
    r = math.isqrt(n)
    max_e = int(math.log(n, 3))

    primes = get_prime_array(n)
    #print(f"Primes({len(primes)}) {primes[0]} to {primes[-1]}")

    # Only interested in p % 8 in (3,5) and odd e
    special_primes = array.array('Q')
    special_primes.extend(filter(lambda p: p % 8 in (3, 5), primes))
    primes = None

    def count_in_ex(n, opi):
        count = n
        if n < special_primes[opi]:
            return count

        for pi in range(opi, len(special_primes)):
            p = special_primes[pi]
            if p > n:
                break

            p2 = p * p
            if p2 > n:
                assert p % 8 in (3, 5)
                count -= n // p
                continue

            # Handle primes with power > 1 or
            pp = p
            for e in range(1, max_e+1, 2):
                assert pp % 8 in (3, 5)
                if pp > n:
                    break

                count -= count_in_ex(n // pp, pi+1)
                # Have to add back all the counts of pp*p
                count += count_in_ex(n // pp // p, pi+1)

                pp *= p2

        return count

    return count_in_ex(n, 0)


def get_three_five_prime_counts(n, primes):
    '''
    Get number of primes % 8 == {3, 5} <= i for important values of i

    Adapted from Lucy_Hedgehog's post in Problem 10
    https://projecteuler.net/thread=10;page=5#111677
    https://math.stackexchange.com/a/2283829/87805
    '''

    r = math.isqrt(n)
    assert primes[-1] >= r

    large_V = [n//i for i in range(1,r+1)]
    small_V = list(range(large_V[-1]-1,0,-1))
    V = large_V + small_V

    # How many numbers <= i of form 8*n + {1, 7} that survive sieving up to and including p
    #Ca = {i: 2 * (i//8) + (i % 8 >= 1) + (i % 8 >= 7)  for i in V}
    # How many numbers <= i of form 8*n + {3, 5} that survive sieving up to and including p
    #Cb = {i: 2 * (i//8) + (i % 8 >= 3) + (i % 8 >= 5)  for i in V}

    C = {
        i:[
            2 * (i//8) + (i % 8 >= 1) + (i % 8 >= 7),
            2 * (i//8) + (i % 8 >= 3) + (i % 8 >= 5),
        ]
        for i in V
    }

    for p in primes:
        if p > r: break
        if p == 2: continue

        c_a, c_b = C[p-1] # count of primes: (8*k + {1,7}, 8*k + {3,5})

        p2 = p*p

        if p % 8 in (1,7):
            for v in V:
                if v < p2: break
                t = C[v // p]
                # count of numbers that don't have divisors less than p
                #   (otherwise all multiplied would already be crossed of)
                #   minus count of primes (which already marked of all multiplies)

                # 1*1 = 1, 1*7 = 7, 7*7 = 1, 7*1 = 7
                #C[v][0] -= Ca[v//p] - c_a
                # 1*3 = 3, 1*5 = 5, 7*3 = 5, 7*5 = 3
                #C[v][1] -= Cb[v//p] - c_b

                C[v][0] -= t[0] - c_a
                C[v][1] -= t[1] - c_b
        else:
            for v in V:
                if v < p2: break
                t = C[v // p]
                # 3*1 = 3, 3*7 = 5, 5*1 = 5, 5*7 = 3
                #Cb[v] -= Ca[v//p] - c_a
                # 3*3 = 1, 3*5 = 7, 5*3 = 7, 5*5 = 1
                #Ca[v] -= Cb[v//p] - c_b

                C[v][1] -= t[0] - c_a
                C[v][0] -= t[1] - c_b

    # Ca also counts 1 which is pseudo "prime"

    '''
    primes = primesieve.primes(n)
    a = [1] + [p for p in primes if p % 8 in (1, 7)]
    b = [p for p in primes if p % 8 in (3, 5)]
    import bisect

    for i in V:
        brute_a = bisect.bisect(a, i)
        brute_b = bisect.bisect(b, i)
        print(i, "\t", Ca[i], Cb[i], "\t", brute_a, brute_b)
        assert Ca[i] == brute_a
        assert Cb[i] == brute_b
    '''
    return {i: t[1] for i, t in C.items()}


def A000047_fast_fast(bits: int) -> int:
    n = 2 ** bits
    r = math.isqrt(n)
    max_e = int(math.log(n, 3))

    # Need slightly more than sqrt(r) primes
    primes = get_prime_array(2 * r + 100)
    #print(f"Primes({len(primes)}) {primes[0]} to {primes[-1]}")

    # Adapted from Lucy_Hedgehog's post in Problem 10
    # https://projecteuler.net/thread=10;page=5#111677
    # https://math.stackexchange.com/a/2283829/87805
    count_special_primes = get_three_five_prime_counts(n, primes)


    # Only interested in p % 8 in (3,5) and odd e
    special_primes = array.array('Q')
    special_primes.extend(filter(lambda p: p % 8 in (3, 5), primes))
    primes = None

    def count_in_ex(n, opi):
        if n < special_primes[opi]:
            return n

        count = n
        for pi in range(opi, len(special_primes)):
            p = special_primes[pi]
            p2 = p * p
            if p2 > n:
                break

            pp = p
            for e in range(1, max_e+1):
                if pp > n:
                    break

                assert pp % 8 in (3, 5)
                count -= count_in_ex(n // pp, pi+1)
                # Have to add back all the counts of pp*p
                count += count_in_ex(n // pp // p, pi+1)

                pp *= p2


        # Handle primes > sqrt(n)
        # for pi in range(pi, len(special_primes)):
        #     p = special_primes[pi]
        #     if p > n:
        #         break
        #     if p % 8 in (3, 5):
        #         count -= n // p

        # Correct for special_primes > pi
        start_p = special_primes[pi]

        next_m = n // start_p - 1
        if next_m == 0:
            break_p = n + 1
        else:
            break_p = n // (next_m + 1) + 1

        for pi in range(pi, len(special_primes)):
            p = special_primes[pi]
            if p >= break_p:
                break
            assert p % 8 in (3, 5)
            count -= n // p

        for m in range(next_m, 0, -1):
            # Count of number of primes with n // p == m
            #   -> Primes in the interval (n // (m + 1), n // m]

            count -= m * (count_special_primes[n // m] - count_special_primes[n // (m+1)])

        return count

    return count_in_ex(n, 0)


def A000047_final(bits: int) -> int:
    n = 2 ** bits
    r = math.isqrt(n)
    max_e = int(math.log(n, 3))

    # Need 2 more primes than sqrt(r) primes
    primes = get_prime_array(r + 1000)
    #print(f"Primes({len(primes)}) {primes[0]} ... {primes[-3:]}")
    assert primes[-2] > r, primes[-5:]

    # Roughly 20-60% of time is taken up with calculating special prime counts
    count_special_primes = get_three_five_prime_counts(n, primes)
    if bits > 35:
        print(f"\tcount_special_primes(2^{bits}) = {count_special_primes[n]}")
    # return count_special_primes[n]

    # Only interested in p % 8 in (3,5) and odd e

    assert primes[-1] < 2 ** 32
    special_primes = array.array('L')
    special_primes.extend(filter(lambda p: p % 8 in (3, 5), primes))
    max_special_prime = special_primes[-1]
    primes = None

    def count_in_ex(n, opi):
        if n < special_primes[opi]:
            return n

        count = n
        max_power = max_e
        for pi in range(opi, len(special_primes)):
            p = special_primes[pi]
            p2 = p * p
            if p2 > n:
                break

            tn = n
            for e in range(1, max_power+1, 2):
                tn //= p
                if tn < p:
                    count -= tn
                    max_power = e
                    break
                count -= count_in_ex(tn, pi+1)

                # Have to add back all the counts of pp*p
                tn //= p
                if tn < p:
                    count += tn
                    max_power = e
                    break
                count += count_in_ex(tn, pi+1)

        # Handle primes > sqrt(n)
        start_p = special_primes[pi]
        assert start_p * start_p > n
        first_m = n // start_p

        last = n // (first_m + 1)
        count_last = count_special_primes[last]
        for m in range(first_m, 0, -1):
            # Count of number of primes with n // p == m
            #   -> Primes in the interval (n // (m + 1), n // m]

            first = last
            last  = n // m

            if m == first_m:
                assert first < last
                assert first <= max_special_prime
                assert first < start_p <= last

            if first < start_p:
                assert m == first_m
                assert start_p <= max_special_prime
                count_first = pi
                #assert count_first = bisect.bisect(special_primes, start_p - 1)
            else:
                # Nice double check of special_prime code
                #count_first = count_special_primes[first]
                count_first = count_last

            count_last = count_special_primes[last]

            assert count_last >= count_first
            count -= m * (count_last - count_first)

        return count

    return count_in_ex(n, 0)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        n = int(sys.argv[1])
        assert n in range(51)

        if n <= 26:
            t0 = time.time()
            count = A000047(n)
            t1 = time.time()
            print(f"A000047          ({n}) = {count}  ({t1-t0:.2f})")

        if n <= 28:
            t0 = time.time()
            count = A000047_fast(n)
            t1 = time.time()
            print(f"A000047_fast     ({n}) = {count}  ({t1-t0:.2f})")

        if n <= 30:
            t0 = time.time()
            count = A000047_fast_fast(n)
            t1 = time.time()
            print(f"A000047_fast_fast({n}) = {count}  ({t1-t0:.2f})")

        t0 = time.time()
        count = A000047_final(n)
        t1 = time.time()
        print(f"A000047_final    ({n}) = {count}  ({t1-t0:.2f})")

    else:
        values = []
        from multiprocessing import Pool
        # Pool too large may run out of memory
        with Pool(6) as p:
            seq = {}
            for n, count in enumerate(p.imap(A000047_final, range(5, 40), chunksize=1), 5):
                seq[n] = count
                print(f"{n}\t{count}")
        print("A000047:", " ".join([str(v) for k, v in sorted(seq.items())]))
