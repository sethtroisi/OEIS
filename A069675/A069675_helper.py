import math
import gmpy2


# How many you want to find
MAX_COUNT = 500

k_count = 3.7               # d = 1000 yields ~264

#for Python
#k_cost = 7 * 1e-10         # d = 1000 takes <XX>s

#for parallel c++
k_cost = 4.14 * 1e-11       # d = 5000 takes ~400s
k_filter_cost = 1.6 * 1e-9  # d = 5000, sieve = 30M takes 10.3s



def optimal_sieve(d, expected_cost):
  all_a_b = d * 81
  non_trivial_a_b = d * 24 # removes 2, 3, 5

  expected_after_sieve = non_trivial_a_b
  cost = 0
  best_cost = expected_cost + 1.0

  prime_pi = 3
  current_prime = gmpy2.mpz(5)
  while True:
    if current_prime < 300000:
      group_size = 1
      if current_prime < 50000:
        current_prime = int(gmpy2.next_prime(current_prime))
      else:
        # when we expect to find the next prime
        current_prime += math.log(current_prime)
    else:
      # do groups of primes at the same time
      group_size = int(current_prime / 100000)
      current_prime += group_size * math.log(current_prime)

    prime_pi += group_size

    expected_after_sieve *= (1 - float(0.98 / current_prime)) ** group_size

    calc_cost = group_size * d * k_filter_cost
    cost += calc_cost

    new_cost = cost + (expected_cost * expected_after_sieve / non_trivial_a_b)

    if new_cost > best_cost:
      break
    best_cost = new_cost

  return (cost,
          expected_cost * expected_after_sieve / non_trivial_a_b,
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
  t_cost = 24 * k_cost * d_cost
  t_count = 24 * k_count * d_count
  return t_cost, t_count


expected_count = 170 # count below a googol
expected_cost = 0

last_print_count = 0

# paired with expected_count = 170 this helps with the initial
# not-quite-so normal zone of the function.
d = 100
while expected_count < MAX_COUNT:
  t_cost, t_count = cost_test_d(d)

  expected_cost += t_cost
  expected_count += t_count

  if int(expected_count) > int(last_print_count):
    #print ("expect {:.0f} around 10^{} (cost: ~~{:.1f} seconds)".format(
    #    expected_count, d, expected_cost))

    sieve_cost, post_sieve_cost, sieve_limit, prime_pi, to_check = \
        optimal_sieve(d, expected_cost)

    print ("expect {:.0f} around 10^{} (optimal sieve: P({}) ~{} leaves {}, cost: ~~{:.2f}) cost: ~~{:.1f} seconds".format(
        expected_count, d, sieve_limit, prime_pi, to_check, sieve_cost, post_sieve_cost))

    last_print_count = expected_count

  d += 1
# '''
