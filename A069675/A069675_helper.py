import math
import gmpy2

# How many you want to find
MAX_COUNT = 450

K_COUNT = 3.7               # d = 1000 yields ~264

#for parallel C++
K_COST = 4.14 * 1e-11       # d = 5000 takes ~400s
K_FILTER_COST = 1.6 * 1e-9  # d = 5000, sieve = 30M takes 10.3s


def optimal_sieve(d, expected_cost):
  non_trivial_a_b = d * 24 # removes 2, 3, 5,

  expected_after_sieve = non_trivial_a_b
  sieve_cost = 0
  best_cost = expected_cost + 1.0

  prime_pi = 3
  current_prime = gmpy2.mpz(5)
  while True:
    if current_prime < 1e5:
      group_size = 1
      current_prime = int(gmpy2.next_prime(current_prime))
    else:
      # do groups of primes at the same time
      group_size = int(current_prime / 100000)
      current_prime += group_size * math.log(current_prime)

    prime_pi += group_size

    filter_rate = (1 - (0.98 / current_prime)) ** group_size
    expected_after_sieve *= filter_rate

    calc_cost = group_size * d * K_FILTER_COST
    sieve_cost += calc_cost

    filter_ratio = expected_after_sieve / non_trivial_a_b
    new_cost = sieve_cost + filter_ratio * expected_cost

    if new_cost > best_cost:
      break
    best_cost = new_cost

  return (sieve_cost,
          expected_cost * filter_ratio,
          int(current_prime),
          prime_pi,
          int(expected_after_sieve))


def cost_test_d(d):
  log_d = d * math.log(10)

  # log_a is trivial compared to log_d
  log_num = log_d # + log_a

  # In theory log_num ^ 2
  # In practice log_num ^ 2.3
  d_cost = log_num ** 2.3
  d_count = 1 / log_num

  # 24 a,b pairs are valid
  t_cost = 24 * K_COST * d_cost
  t_count = 24 * K_COUNT * d_count
  return t_cost, t_count

def maybe_M(n):
    if n < 1e7:
        return n
    if n < 1e9:
        return "{:.1f}M".format(n / 1e6)
    if n < 1e12:
        return "{:.1f}B".format(n / 1e9)
    return "{:.1f}T".format(n / 1e12)

def maybe_H(n):
    if n < 3 * 3600:
        return "{:.1f} seconds".format(n)
    if n < 2 * 86400:
        return "{:.1f} hours".format(n / 3600.0)
    if n < 365 * 86400:
        return "{:.1f} days".format(n / 86400.0)
    return "{:.1f} years".format(n / 86400.0 / 365.0)

expected_count = 170 # count below a googol
expected_cost = 0

last_print_count = 0

# paired with expected_count = 170 this helps with the initial
# not-quite-so normal zone of the function.
d = 100
while expected_count < MAX_COUNT:
  mult = 1 if d < 1000 else 10
  t_cost, t_count = cost_test_d(d)

  expected_cost += mult * t_cost
  expected_count += mult * t_count

  if int(expected_count) > int(last_print_count):
    sieve_cost, post_sieve_cost, sieve_limit, prime_pi, to_check = \
        optimal_sieve(d, expected_cost)

    sieve_stats = "optimal sieve: PrimePi({}) ~= {}, leaves {} cost ~~{}".format(
        maybe_M(sieve_limit), maybe_M(prime_pi),
        to_check,
        maybe_H(sieve_cost))

    print ("expect {:.0f} around 10^{} ({}) cost: ~~{}".format(
        expected_count, d, sieve_stats, maybe_H(post_sieve_cost)))

    last_print_count = expected_count

  d += mult

