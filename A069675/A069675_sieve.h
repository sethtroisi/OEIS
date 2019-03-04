#include <atomic>

#define START_DIGIT 1
#define MAX_DIGITS  200000

#define ONE_MILLION 1000000L
#define ONE_BILLION 1000000000L
#define SIEVE_LIMIT (10000L * ONE_MILLION)

extern long is_prime[MAX_DIGITS+1][10][10];
