import math

# ignore small numbers where primes not as uniformaly distributed

# Total primes <= 10^21
total = 13064499

# Count of primes 10^19..10^20, 10^20..10^21
# these are 20 digit and 21 digit primes
count = [3492301, 3332745]
start_digit = 22

for digits in range(start_digit, 110):
    # Building numbers of "digits" digits

    MIN = 10 ** (digits - 1)
    MAX = 10 ** (digits)

    # intervals of last generation, assumes all numbers are equal likely to be
    # prime (this is an okay approximation)

    intervals = []
    for interval in range(1, 1000):
      t = interval / 1000 * (MAX / 100)
      intervals.append(t)

    count_last = 0

    old_count = count.pop(0)
    new_count = 0
    for left in range(1, 10):
        right_count = 10

        for interval in intervals:
          rough_count = old_count / len(intervals)
          prob_has_child = right_count / math.log(left * MIN + interval)

          new_count += rough_count * prob_has_child

    new_count = int(new_count)
    total += new_count

    print ("d: {}\tnew: {:<9d} (old: {:<9d}, r: {:.2f})\ttotal: {}".format(
      digits, new_count, old_count, new_count / (old_count + 1), total))

    count.append(new_count)
