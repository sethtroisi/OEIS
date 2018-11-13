#define START_DIGIT 1
#define MAX_DIGITS  3000
#define MAX_DIGITS_P1 (MAX_DIGITS + 1)

#define ONE_MILLION 1000000L
#define SIEVE_LIMIT 2 * ONE_MILLION

//#define PRIME_PI_1M 78498
#define PRIME_PI_1M 70000

#ifdef __CUDACC__
    #define CUDAGLOBAL __global__
    #define CUDAHOSTDEV __host__ __device__
#else
    #define CUDAGLOBAL
    #define CUDAHOSTDEV
#endif

CUDAHOSTDEV bool test_p(
//    long is_prime_ref[MAX_DIGITS][10][10],
    long *is_prime_ref,
    long p,
    long mods[24]);
