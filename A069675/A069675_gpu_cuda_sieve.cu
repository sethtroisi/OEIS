#include <iostream>

#include "A069675_gpu_shared.h"

using namespace std;

#define checkCuda(expr) {                          \
  auto status = (expr);                             \
  if (status != cudaSuccess) {                       \
    cerr << "cuda Error on line " << __LINE__ << ": " \
         << cudaGetErrorString(status) << endl;        \
    exit(EXIT_FAILURE);                                 \
  }                                                      \
}


// TODO figure out what to set here
#define GRID_SIZE 16
#define BLOCK_SIZE 4

CUDAGLOBAL void FilterSieveKernel(
    long *is_prime_ref,
    long *div_mods,
    long *primes,
    long pi_start,
    long pi_end,
    bool *results) {

  // TODO LOTS of cuda calls
  int index = threadIdx.x + blockIdx.x * BLOCK_SIZE;

  long range = pi_end - pi_start;
  long bottom = range * index / (BLOCK_SIZE * GRID_SIZE);
  long top = range * (index + 1) / (BLOCK_SIZE * GRID_SIZE);
  bottom += pi_start;
  top += pi_start;

  for (long pi = bottom; pi < top; pi++) {
    results[pi] = test_p(is_prime_ref, primes[pi], div_mods + 24 * pi);
  }
}

void FilterSieveKernelHost(
    long is_prime_ref[MAX_DIGITS_P1][10][10],
    long div_mods[][24],
    long primes[],
    long pi_start,
    long pi_end,
    bool results[]) {

  // Probably beneficial to compress 100 bits vs 6400
  long *d_is_prime = NULL;
  checkCuda(cudaMalloc(&d_is_prime, sizeof(long) * MAX_DIGITS_P1 * 10 * 10));
  checkCuda(cudaMemcpy(d_is_prime, (long*)is_prime_ref, sizeof(long) * MAX_DIGITS_P1 * 10 * 10, cudaMemcpyHostToDevice));

  long *d_div_mods = NULL;
  checkCuda(cudaMalloc(&d_div_mods, sizeof(long) * 24 * pi_end));
  checkCuda(cudaMemcpy(d_div_mods, div_mods, sizeof(long) * 24 * pi_end, cudaMemcpyHostToDevice));

  long *d_primes = NULL;
  checkCuda(cudaMalloc(&d_primes, sizeof(long) * pi_end));
  checkCuda(cudaMemcpy(d_primes, primes, sizeof(long) * pi_end, cudaMemcpyHostToDevice));

  bool *d_results = NULL;
  checkCuda(cudaMalloc(&d_results, pi_end));
  checkCuda(cudaMemset(d_results, 0, pi_end));

  FilterSieveKernel<<<GRID_SIZE, BLOCK_SIZE>>>(
      d_is_prime, d_div_mods, d_primes, pi_start, pi_end, d_results);

  checkCuda(cudaMemcpy(results, d_results, pi_end, cudaMemcpyDeviceToHost));

  cudaFree(d_is_prime);
  cudaFree(d_div_mods);
  cudaFree(d_primes);
  cudaFree(d_results);
}

