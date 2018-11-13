
#include <cuda.h>
#include <cuda_runtime.h>

#include "A069675_gpu_shared.h"

// TODO figure out what to set here
#define GRID_SIZE 32
#define BLOCK_SIZE 32

CUDAGLOBAL void FilterSieveKernel(
    long is_prime_ref[MAX_DIGITS][10][10],
    long div_mods[][24],
    long primes[],
    long pi_start,
    long pi_end,
    bool results[]);

