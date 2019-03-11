// Include Guard
#pragma once

#define START_DIGIT 1
#define MAX_DIGITS  1000000

#define ONE_MILLION 1000000L
#define ONE_BILLION 1000000000L
//#define SIEVE_LIMIT (1000000L * ONE_MILLION)
#define SIEVE_LIMIT (100L * ONE_MILLION)

// p = divisibile by p
// 0 = (to test)
// -1 = Prime
// -2 = Unknown Composite
long is_prime[MAX_DIGITS+1][10][10] = {};
