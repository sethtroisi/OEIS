"""
Clear small composites from (3^k - 7) / 2
"""

import array
import math

import primesieve
import sympy


def bounded_discrete_log(prime, a, b, BOUND):
    """
    Same as sympy.ntheory.residue_ntheory.discrete_log

    Stops (and returns None) if discrete logarithm would be greater than BOUND
    Uses Baby-step Giant-step algorithm
    """

    # Requires factoring prime-1 but shrug
    #BOUND = min(b_order, BOUND)
    BOUND = min(prime-1, BOUND)

    sqrt = math.isqrt(BOUND) + 1
    T = {}

    power = 1
    for i in range(sqrt):
        T[power] = i
        power = (power * b) % prime

    backwards_sqrt_steps = pow(b, -sqrt, prime)
    assert (pow(b, sqrt, prime) * backwards_sqrt_steps) % prime == 1

    power = a
    for i in range(sqrt):
        if power in T:
            return i * sqrt + T[power]

        # skip "backwards" m steps
        power = (power * backwards_sqrt_steps) % prime

    return None


def sieve_terms(status, prime_limit):
    """
    Handle 2, 3, 5, 7 so that loop starts on

    2 divides  <=>  k % 2 == 1
    3 never divides
    5 divides  <=>  k % 4 == 3 (handled by 2)
    7 never divides
    """

    for k in range(3, len(status), 2):
        assert pow(3, k, 4) == (7 % 4)
        status[k] = 0

    for prime in primesieve.primes(7+1, prime_limit):
        # Find first prime that this divides
        #    <=>
        # 3^k = 7 MOD (2*prime)

        '''
        goal = 7 % prime
        power_three = 3
        for i in range(1, prime):
            if power_three == goal:
                assert (pow(3, i, prime) - 7) % prime == 0
                for j in range(i, N+1, prime-1):
                    assert (pow(3, j, prime) - 7) % prime == 0
                    status[j] = 0
                break

            power_three = (power_three * 3) % prime
        '''
        '''
        try:
            first = sympy.ntheory.residue_ntheory.discrete_log(prime, 7, 3)
            if first <= len(status):
                assert first == bounded_discrete_log(prime, 7, 3, min(prime, len(status)))
            assert pow(3, first, prime) == 7
            assert pow(3, first + (prime-1), prime) == 7
            for j in range(first, len(status), prime-1):
                #assert (pow(3, j, prime) == 7
                status[j] = 0
        except ValueError as e:
            # Log does not exist OR two numbers should be relatively prime
            #print("Failed for", prime, "\t", e)
            pass
        except KeyboardInterrupt:
            print(f"Stopping during prime={prime}")
            return
        '''

        b_order = sympy.ntheory.residue_ntheory.n_order(3, prime)
        BOUND = min(b_order, len(status))
        first = bounded_discrete_log(prime, 7, 3, BOUND)
        if first:
            #z = sympy.ntheory.residue_ntheory.discrete_log(prime, 7, 3)
            #assert first == z, (prime, first, z)
            assert pow(3, first, prime) == 7
            assert pow(3, b_order, prime) == 1
            assert pow(3, first + b_order, prime) == 7

            for j in range(first, len(status), b_order):
                #assert (pow(3, j, prime) == 7
                status[j] = 0

        if prime % 100000 <= 3:
            print(f"\t{prime} {sum(status)} remaining")



def main():
    N = 30000
    PRIME_LIMIT = 10 ** 6
    print(f"Sieving up to (3^k-7)/2 k=2..{N} with primes <= {PRIME_LIMIT:,}")
    print()

    status = array.array('B', [1]) * (N+1)
    status[0] = status[1] = 0

    sieve_terms(status, PRIME_LIMIT)

    print("Writing", sum(status), "terms to disk")
    with open("to_test.txt", "w") as f:
        for i, s in enumerate(status):
            if s:
                f.write(f"(3^{i}-7)/2\n")

if __name__ == "__main__":
    main()
