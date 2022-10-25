import array
import math
import sys
import tqdm

def gen_primes(N):
  is_prime = array.array('B', [1]) * (N+1)
  is_prime[0] = is_prime[1] = 0
  for p in range(3, math.isqrt(N)+1, 2):
    if is_prime[p]:
      for m in range(p*p, N+1, 2*p):
        is_prime[m] = 0

  primes = [2]
  for p in range(3, N+1, 2):
    if is_prime[p]:
      primes.append(p)

  return primes

def tonelli_shanks(a, p):
  q = p - 1
  s = 0
  while q&1 == 0: # q % 2 == 0
    q >>= 1
    s += 1

  if s == 1:  # n % 4 == 3
    # Easy euler criterion case
    assert p % 4 == 3
    return pow(a, (p+1)//4, p)

  # Find a non-residual
  for z in range(2, p):
    if pow(z, (p-1)//2, p) == (p-1):
      break

  b2 = pow(z, q, p)
  r = pow(a, (q+1)//2, p)
  t = pow(a, q, p)
  m = s
  t2 = 0
  # Invariant
  while (t - 1) % p != 0:
    # Find the first power that's a divisor
    t2 = t
    for i in range(m):
      if (t2 - 1) % p == 0:
        break
      t2 = (t2 * t2) % p

    b = pow(b2, 2 ** (m - i - 1), p)
    r = (r * b) % p
    b2 = (b * b) % p
    t = (t * b2) % p
    m = i
  return r


def factor_offset(N, twice_square, offset, primes):
  """
  Factor a^2 - offset  <=>  a^2 = offset mod p

  offset is negative so that it's the quadratic residual

  offset is a small odd number and can be negative
  """
  assert (offset == 0) or (offset % 2 == 1)
  square_mult = 2 if twice_square else 1

  # factored, sigma
  status = array.array('L', [1, 1]) * (N+1)
  factor_count = 0

  def remove_factor(base, prime):
    nonlocal factor_count

    for index in range(base, N+1, prime):
      num = square_mult * index * index - offset
      num, m = divmod(num, prime)
      assert m == 0

      pp = prime
      exponent = 1

      # 20-40% of time is spent checking for p^2
      while num:
        num, m = divmod(num, prime)
        if m != 0:
          break
        pp *= prime
        exponent += 1

      factor_count += 1  # exponent?
      #assert (index * index - offset) % pp == 0
      #print("\t\t", prime, index)
      status[2*index+0] *= pp
      status[2*index+1] *= (prime * pp - 1) // (prime - 1)

  start_index = min(i for i in range(100) if square_mult * i * i - offset > 2)

  for prime in primes[1:]:
    if twice_square:
      # 2*a*a = offset
      # a*a = offset/2
      inv_two = (prime+1)//2
      #assert (2 * inv_two) % p == 1
      residual = (inv_two * offset) % prime
    else:
      residual = offset % prime

    if residual == 0:
        # first multiple of prime >= start_index
        start = ((start_index + prime - 1) // prime) * prime
        remove_factor(start, prime)
    else:
      # Check if quadratic residue exists
      legendre = pow(residual, (prime-1)//2, prime)

      if legendre == 1:
        # Find "square root" of offset, using tonelli-shanks
        base = tonelli_shanks(residual, prime)
        assert 1 <= base < prime
        #assert (square_mult * base * base - offset) % prime == 0

        # Want to find the first time when divisibe by prime^2
        # See: https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm#Tonelli's_algorithm_will_work_on_mod_p^k
        # c = offset = residual
        # x = base
        #prime_2 = prime * prime
        #test = pow(base, prime, prime_2) * pow(residual, (prime_2 - 2 * prime + 1) // 2, prime_2) % prime_2
        #assert (square_mult * test * test - residual) % prime_2 == 0

        assert base != prime - base
        for b in [base, prime - base]:
          # sometimes need to adjust to first positive index
          if b < start_index:
            b += prime

          remove_factor(b, prime)

  print(f"\t{factor_count:,} prime factors in {square_mult} * n^2 {-offset:+}")
  #print()
  #print("results")
  for i in range(start_index, N+1):
    num = square_mult * i * i - offset
    rem = num // status[2*i+0]

    # Handle 2's
    if (rem & 1) == 0:
      twos = 0
      while rem & 1 == 0:
        rem >>= 1
        twos += 1

      status[2*i+1] *= (1 << (twos + 1)) - 1

    if rem > 1:
      assert rem > N and rem % 2 == 1, (i, num, rem)
      status[2*i+1] *= (1 + rem)

    #print("\t", i, num, "\t", rem, "\t", status[2*i+1])

  #print(f"\tSigmas for {square_mult}*n^2 {-offset:+} computed", status[1:10:2])
  return status[1::2]


if __name__ == "__main__":
  N = 10 ** 12 if len(sys.argv) != 2 else int(sys.argv[1])
  print(f"Testing up to {N:,} ({N:.1e})")

  STOP = math.isqrt(N)
  HALF_STOP = math.isqrt(N//2)
  primes = gen_primes(STOP)
  half_primes = [p for p in primes if p <= HALF_STOP]

  print(f"primes({STOP:,}) = {len(primes)} {primes[:5]} ... {primes[-2:]}")

  for twice_square in (False, True):
    def factor_offset_bounded(offset):
      return factor_offset(
        HALF_STOP if twice_square else STOP,
        twice_square,
        offset,
        half_primes if twice_square else primes)

    sigmas_squares = factor_offset_bounded(offset=0)

    OFFSETS = [1,-1, 3, -3, 5, -5, 7, -7,9, -9, 11, -11, 13, -13, 15, -15, 17, -17, 19, -19]
    for offset in OFFSETS:
      sigmas_offset = factor_offset_bounded(offset)

      square_mult = (1 + twice_square)
      start_i = min(i for i in range(2, 100) if square_mult*i*i-offset > 2)
      for i, sigma_a in enumerate(sigmas_squares[start_i:], start_i):
          num = square_mult * i * i
          sigma_b = sigmas_offset[i]
          if sigma_a - offset == sigma_b:
            print()
            print(f"MATCH @ {i} = {num} and {num - offset} | {sigma_a} {-offset:+} vs {sigma_b}")
            print()
