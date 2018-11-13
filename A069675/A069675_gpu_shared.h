#define START_DIGIT 1
#define MAX_DIGITS  5000

#define ONE_MILLION 1000000L
#define SIEVE_LIMIT 1 * ONE_MILLION

// TODO figure out what to set here
#define GRID_SIZE 32
#define BLOCK_SIZE 32

#ifdef __CUDACC__
    #define CUDAGLOBAL __global__
    #define CUDAHOSTDEV __host__ __device__
#else
    #define CUDAGLOBAL
    #define CUDAHOSTDEV
#endif

CUDAHOSTDEV bool test_p(
    long is_prime_ref[MAX_DIGITS][10][10],
    long p,
    long mods[24]);
