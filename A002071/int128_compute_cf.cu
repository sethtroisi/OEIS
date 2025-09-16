/* int128_compute_cf.h: header for (GPU) based continued fractions.

  Copyright 2025 Seth Troisi

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
more details.

You should have received a copy of the GNU General Public License
along with this program; see the file COPYING.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#ifndef __CUDACC__
#error "This file should only be compiled with nvcc"
#endif

#include "int128_compute_cf.h"

#include <cassert>
#include <cstdio>


inline void cuda_check(cudaError_t status, const char *action=NULL, const char *file=NULL, int32_t line=0) {
  if (status != cudaSuccess) {
    fprintf (stderr, "CUDA error (%d) occurred: %s\n", status, cudaGetErrorString(status));
    if (action!=NULL)
      fprintf (stderr, "While running %s   (file %s, line %d)\n", action, file, line);
    exit(EXIT_FAILURE);
  }
}

#define CUDA_CHECK(action) cuda_check(action, #action, __FILE__, __LINE__)

// kernel implementation using __int128 data types
__global__ void kernel_iterate_cf(
        uint64_t MAX_CF,
        uint32_t count,
        __uint128_t *data,
        uint32_t *results
) {
  // decode an instance_i number from the blockIdx and threadIdx
  int32_t instance_i = blockIdx.x*blockDim.x + threadIdx.x;
  if(instance_i >= count)
    return;

  __uint128_t x, a0, two_a0, a, b, c;
  __uint128_t t;

  { // Setup
    x  = data[6*instance_i+0];
    a0 = data[6*instance_i+1];
    two_a0 = data[6*instance_i+2];
    b  = data[6*instance_i+3];
    c  = data[6*instance_i+4];
    a  = data[6*instance_i+5];
  }

  /*
  if (instance_i == 0) {
    printf("GPU %2d | setup %lu %lu %lu | %lu %lu %lu\n", instance_i,
            (uint64_t) x, (uint64_t) a0, (uint64_t) two_a0,
            (uint64_t) b, (uint64_t) c, (uint64_t) a);
  }
  */
  results[instance_i] = 0;

  for (uint64_t i = 2; i <= MAX_CF; i++) {
    // b = a*c - b;
    t = a * c;
    b = t - b;

    // c = (x - b*b) / c
    t = b * b;
    t = x - t;

    if (c == 1) {
        c = t;
        a = a0 + b;
    } else {
        c = t / c;
        t = a0 + b;
        a = t / c;
    }

    /*
    if (instance_i == 0) {
      printf("GPU %2d @ %3lu | %lu -> goal: %lu | %lu %lu %lu\n", instance_i, i,
              (uint64_t) x, (uint64_t) two_a0, (uint64_t) b, (uint64_t) c, (uint64_t) a);
    }
    */

    // Can this be moved under the c == 1 branch?
    // done = (a == two_a0)
    if ((two_a0 == a) && results[instance_i] == 0) {
      results[instance_i] = i + 1;
    }
  }
}



void PessemisticCf::store_data(vector<pair<__uint128_t, __uint128_t>>& D_a0) {
    size_t    count = D_a0.size();
    assert( count <= max_D_count );
    for (size_t i = 0; i < count; i++) {
        auto [D, a0] = D_a0[i];
        auto c = D - a0 * a0;

        data[6*i + 0] = D;
        data[6*i + 1] = a0;
        data[6*i + 2] = a0 << 1; // two_a0
        data[6*i + 3] = a0; // b
        data[6*i + 4] = c;  // c
        data[6*i + 5] = (a0 << 1) / c; // a
        //printf("%2lu | %lu %lu %lu\n", i, (uint64_t) D, (uint64_t) a0, (uint64_t) (a0 << 1));
    }
}


PessemisticCf::PessemisticCf(size_t D_size) {
    CUDA_CHECK(cudaEventCreate (&start));
    CUDA_CHECK(cudaEventCreate (&stop));

    max_D_count = D_size;
    data_size = 6 * sizeof(__uint128_t) * D_size;
    data = (__uint128_t*) malloc(data_size);

    CUDA_CHECK(cudaMalloc((void **)&gpu_data, data_size));
    CUDA_CHECK(cudaMalloc((void **)&gpu_results, sizeof(uint32_t) * D_size));
}

PessemisticCf::~PessemisticCf() {
    free(data);
    CUDA_CHECK(cudaFree(gpu_data));
    CUDA_CHECK(cudaFree(gpu_results));
    CUDA_CHECK(cudaEventDestroy (start));
    CUDA_CHECK(cudaEventDestroy (stop));
}

void PessemisticCf::run(size_t MAX_CF, vector<pair<__uint128_t, __uint128_t>>& D_a0, vector<uint32_t> &valid, int verbose) {
    size_t    count = D_a0.size();

    int32_t   IPB = 128;       // instances per block
    size_t    BLOCK_COUNT;     // How many blocks to cover all count
    BLOCK_COUNT = (count + IPB - 1) / IPB;

    //kernel_info((const void*)kernel_iterate_cf, true);

    store_data(D_a0);

    CUDA_CHECK(cudaEventRecord (start));

    // Copy data
    if (verbose)
      printf("Copying %'lu bytes of data to GPU\n", data_size);
    CUDA_CHECK(cudaMemcpy(gpu_data, data, data_size, cudaMemcpyHostToDevice));

    if (verbose)
      printf("%lu -> GPU_INT128 running kernel<%lu block x %d threads>\n",
          count, BLOCK_COUNT, IPB);

    /* Call CUDA Kernel. */
    kernel_iterate_cf<<<BLOCK_COUNT, IPB>>>(MAX_CF, count, gpu_data, gpu_results);

    /* sync the device */
    CUDA_CHECK(cudaDeviceSynchronize());

    float gputime;
    CUDA_CHECK(cudaEventRecord (stop));
    CUDA_CHECK(cudaEventSynchronize (stop));
    cudaEventElapsedTime (&gputime, start, stop);

    // Copy data back from GPU memory
    if (verbose)
      printf("Copying results back to CPU ...\n");
    CUDA_CHECK(cudaMemcpy(static_cast<void*>(valid.data()), gpu_results, sizeof(uint32_t) * count, cudaMemcpyDeviceToHost));
}

#ifdef __CUDA_ARCH__
  #if __CUDA_ARCH__ < 350
    #error "Unsupported architecture"
  #endif
#endif
