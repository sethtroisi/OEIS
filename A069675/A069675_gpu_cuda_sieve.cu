#include <array>
#include <cassert>
#include <cstdio>
#include <vector>

#include "A069675_gpu_shared.h"

using namespace std;

// TODO figure out what to set here
#define GRID_SIZE 32
#define BLOCK_SIZE 32

CUDAGLOBAL void FilterSieveKernel(
    long is_prime_ref[MAX_DIGITS][10][10],
    long div_mods[][24],
    long primes[],
    long pi_start,
    long pi_end,
    bool results[]) {

  //int index = threadIdx.x + blockIdx.x * BLOCK_SIZE;

  for (long pi = pi_start; pi < pi_end; pi++) {
    test_p(is_prime_ref, primes[pi], div_mods[pi]);
  }
}

