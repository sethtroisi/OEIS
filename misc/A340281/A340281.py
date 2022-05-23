import math
import sys

import gmpy2
import primesieve

# See Seth Troisi's ProjectEuler countsig code.

def gen_all_signatures_with_sigma(sigma, min_d=1):
    # Generate all signatures such that (e1+1) * (e2+1) * (e3+1) ... = sigma

    assert sigma > 1
    assert sigma < 10000, "Replace this code with sympy.divisors()"


    divisors = [sigma]
    for d in range(2, math.isqrt(sigma)+1):
        if sigma % d == 0:
            divisors.append(d)
            divisors.append(sigma // d)

    # choose e1
    for e1 in sorted(divisors):
        if e1 < min_d:
            continue

        if e1 == sigma:
            yield (e1-1,)
        else:
            for sig in gen_all_signatures_with_sigma(sigma // e1, min_d=e1):
                yield (e1-1,) + sig


def nth_root(N, nth):
    if nth == 2:
        return math.isqrt(N)
    root = max(1, int(N**(1 / nth)) - 100)
    assert root ** nth <= N, (N, root, nth, root ** nth)
    while (root + 1) ** nth <= N:
        root += 1
    return root


def enumerate_signature_A340281(signature, N):
    """Enumerate all odd numbers with sig=signature <= N."""

    # Odd divisors -> no factors of 2

    # TODO determine max value needed here sqrt(N)?
    primes = primesieve.primes(3, 10**6)

    used_primes=set()
    def enumerate_sig_slow(n, sig_index, min_prime_index):
        exp = signature[sig_index]
        allow_new = (sig_index == 0) or signature[sig_index-1] > exp
        min_prime_index = 0 if allow_new else min_prime_index

        if sig_index == len(signature) - 1:
            for pi in range(min_prime_index, len(primes)):
                prime = primes[pi]
                if prime in used_primes:
                    continue
                prime_power = prime ** signature[sig_index]
                if prime_power > n:
                    break

                yield prime ** signature[sig_index]
            return

        for pi in range(min_prime_index, len(primes)):
            prime = primes[pi]
            prime_power = prime ** signature[sig_index]
            temp = n // prime ** signature[sig_index]
            if temp == 0:
                break


            if prime not in used_primes:
                used_primes.add(prime)
                for z in enumerate_sig_slow(temp, sig_index + 1, pi + 1):
                    yield prime_power * z
                used_primes.remove(prime)

    for z in enumerate_sig_slow(N, 0, 0):
        yield z


def smallest_signature_with_next_prime(signature, found):
    assert signature == tuple(sorted(signature, reverse=True)), signature

    def find_smaller(n):
        """Find a smaller n that is prime and matches signature"""
        for sig_n in enumerate_signature_A340281(signature, n):
            assert sig_n % 2 == 1
            # search for is_prime(2^e * sig_n + 1) less than found
            test = sig_n
            while True:
                test *= 2
                if test >= n:
                    break

                if gmpy2.is_prime(test + 1):
                    two_power = int(math.log2(test // sig_n))
                    print(f"\t\tis_prime(2^{two_power} * {sig_n} + 1 = {test}) = True")
                    assert test < n
                    return test + 1

    # Iterate signature looking for n+1 = prime
    while True:
        print(f"\t\tSearching for Num[{signature}] <= {found}")
        smaller = find_smaller(found-1)
        if smaller:
            assert smaller < found
            found = smaller
        else:
            # Didn't find any new best
            break

    return found


def A340281(N):
    """Build all terms up to N."""

    if N <= 0:
        assert False, N

    if N <= 2:
        return N + 1

    LIMIT = 10 ** 200
    smallest = LIMIT
    for si, signature in enumerate(gen_all_signatures_with_sigma(N-1)):
        smallest_with_sig = smallest_signature_with_next_prime(signature[::-1], smallest)
        smallest = min(smallest, smallest_with_sig)
        print(f"\t{si:<4} signature: {signature!s:30} -> {smallest_with_sig}/{smallest}")

    if smallest == LIMIT:
        print("FAILED TO FIND EXAMPLE")
        exit(1)

    return smallest



if __name__ == "__main__":
    with open("b340281.txt", "w") as f:
        for n in range(1, 200+1):
            an = A340281(n)
            f.write(f"{n} {an}\n")
            print(n, an)
