import array
import bisect
import math
import primesieve

def A000047(n: int) -> int:
    # range = [0, last), manually correct for 2 ** n at the end
    last  = 2 ** n

    primes = primesieve.primes(3, last)
    #print(f"Primes {primes[0]} to {primes[-1]}")

    # easier than array
    status = bytearray(last)
    for i in range(last):
        status[i] = 1

    max_e = n
    for p in primes:
        pp = 1
        for e in range(1, max_e+1):
            pp *= p
            if pp > last:
                max_e = e
                break

            if pp % 8 in (3, 5):
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

    # -1 for 0, +1 for 2^n
    count = sum(status)
    #for i, s in enumerate(status):
    #    if s: print(i)
    return count


def A000047_fast(bits: int) -> int:

    n = 2 ** bits
    r = math.isqrt(n)
    max_e = int(math.log(n, 3))

    # Need slightly more than sqrt(r) primes
    primes = primesieve.primes(n) #3 * r)
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
    r = math.isqrt(n)
    assert primes[-1] >= r

    V = [n//i for i in range(1,r+1)]
    V += list(range(V[-1]-1,0,-1))

    # How many numbers <= i of form 8*n + {1, 7} that survive sieving up to and including p
    Ca = {i: 2 * (i//8) + (i % 8 >= 1) + (i % 8 >= 7)  for i in V}
    # How many numbers <= i of form 8*n + {3, 5} that survive sieving up to and including p
    Cb = {i: 2 * (i//8) + (i % 8 >= 3) + (i % 8 >= 5)  for i in V}

    for p in primes:
        if p > r: break
        if p == 2: continue

        c_a = Ca[p-1] # count of primes, 8*k + {1,7}
        c_b = Cb[p-1] # count of primes, 8*k + {3,5}

        p2 = p*p
        p3 = p2*p
        for v in V:
            if v < p2: break
            if p % 8 in (1,7):
                # count of numbers that don't have divisors less than p
                #   (otherwise all multiplied would already be crossed of)
                #   minus count of primes (which already marked of all multiplies)

                # 1*1 = 1, 1*7 = 7, 7*7 = 1, 7*1 = 7
                Ca[v] -= Ca[v//p] - c_a

                # 1*3 = 3, 1*5 = 5, 7*3 = 5, 7*5 = 3
                Cb[v] -= Cb[v//p] - c_b
            else:
                # 3*1 = 3, 3*7 = 5, 5*1 = 5, 5*7 = 3
                Cb[v] -= Ca[v//p] - c_a

                # 3*3 = 1, 3*5 = 7, 5*3 = 7, 5*5 = 1
                Ca[v] -= Cb[v//p] - c_b

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
    return Cb


def A000047_fast_fast(bits: int) -> int:
    n = 2 ** bits
    r = math.isqrt(n)
    max_e = int(math.log(n, 3))

    # Need slightly more than sqrt(r) primes
    primes = primesieve.primes(3 * r)
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

    # Need slightly more than sqrt(r) primes
    primes = primesieve.primes(3 * r)
    #print(f"Primes({len(primes)}) {primes[0]} to {primes[-1]}")

    # Adapted from Lucy_Hedgehog's post in Problem 10
    # https://projecteuler.net/thread=10;page=5#111677
    # https://math.stackexchange.com/a/2283829/87805
    count_special_primes = get_three_five_prime_counts(n, primes)

    special_primes = array.array('Q')
    special_primes.extend(filter(lambda p: p % 8 in (3, 5), primes))
        # Only interested in p % 8 in (3,5) and odd e
    max_special_prime = special_primes[-1]
    primes = None

    def count_in_ex(n, opi):
        if n < special_primes[opi]:
            return n

        count = n
        last_processed = 0
        max_power = max_e
        # Only interested in p % 8 in (3,5) and odd e
        for pi in range(opi, len(special_primes)):
            p = special_primes[pi]
            p2 = p * p
            if p2 > n:
                break

            pp = 1
            for e in range(1, max_power+1):
                pp *= p
                if pp > n:
                    max_power = e
                    break

                if pp % 8 in (3, 5):
                    last_processed = p
        #            print("\t", n, opi, " ", count, "\t", p, e, " ", pp, "\t", n // pp)
                    count -= count_in_ex(n // pp, pi+1)
                    # Have to add back all the counts of pp*p
                    count += count_in_ex(n // pp // p, pi+1)

        # Handle primes > sqrt(n)
        start_p = special_primes[pi]
        assert start_p * start_p > n
        first_m = n // start_p

        for m in range(first_m, 0, -1):
            # Count of number of primes with n // p == m
            #   -> Primes in the interval (n // (m + 1), n // m]

            first = n // (m + 1)
            last  = n // m

            if m == first_m:
                assert first < last
                assert first <= max_special_prime
                assert first < start_p <= last

            if first < start_p:
                assert m == first_m
                assert start_p <= max_special_prime
                count_first = bisect.bisect(special_primes, start_p - 1)
            else:
                # Nice double check of special_prime code
                count_first = count_special_primes[first]
                #test = bisect.bisect(special_primes, first)
                #assert count_first == test, (count_first, test)

            count_last = count_special_primes[last]
            #test = bisect.bisect(special_primes, last)
            #assert count_last == test

            assert count_last >= count_first

            count -= m * (count_last - count_first)

        return count

    return count_in_ex(n, 0)


if __name__ == "__main__":
    n = 35
    if n <= 25:
        count = A000047(n);
        print(f"A000047          ({n}) = {count}")

    if n <= 28:
        count = A000047_fast(n);
        print(f"A000047_fast     ({n}) = {count}")
        print()

    if n <= 30:
        count = A000047_fast_fast(n);
        print(f"A000047_fast_fast({n}) = {count}")

    count = A000047_final(n);
    print(f"A000047_final    ({n}) = {count}")

