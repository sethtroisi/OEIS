/* cgbn_compute_cf.h: header for CGBN (GPU) based continued fractions.

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

#include "cgbn_compute_cf.h"

#include <cassert>

// GMP import must proceed cgbn.h
#include <gmp.h>
#include <cgbn.h>
#include <cuda.h>

// Can dramatically change compile time
#if 1
    #define FORCE_INLINE __forceinline__
#else
    #define FORCE_INLINE
#endif

#define CHECK_ERROR true


inline void cuda_check(cudaError_t status, const char *action=NULL, const char *file=NULL, int32_t line=0) {
  if (status != cudaSuccess) {
    fprintf (stderr, "CUDA error (%d) occurred: %s\n", status, cudaGetErrorString(status));
    if (action!=NULL)
      fprintf (stderr, "While running %s   (file %s, line %d)\n", action, file, line);
    exit(EXIT_FAILURE);
  }
}

#define CUDA_CHECK(action) cuda_check(action, #action, __FILE__, __LINE__)

// support routine copied from  "CGBN/samples/utility/support.h"
void cgbn_check(cgbn_error_report_t *report, const char *file=NULL, int32_t line=0) {
  // check for cgbn errors

  if(cgbn_error_report_check(report)) {
    fprintf (stderr, "\n");
    fprintf (stderr, "CGBN error occurred: %s\n", cgbn_error_string(report));

    if(report->_instance!=0xFFFFFFFF) {
      fprintf (stderr, "Error reported by instance %d", report->_instance);
      if(report->_blockIdx.x!=0xFFFFFFFF)
        fprintf (stderr, ", blockIdx=(%d, %d, %d)", report->_blockIdx.x, report->_blockIdx.y, report->_blockIdx.z);
      if(report->_threadIdx.x!=0xFFFFFFFF)
        fprintf (stderr, ", threadIdx=(%d, %d, %d)", report->_threadIdx.x, report->_threadIdx.y, report->_threadIdx.z);
      fprintf (stderr, "\n");
    }
    else {
      fprintf (stderr, "Error reported by blockIdx=(%d %d %d)", report->_blockIdx.x, report->_blockIdx.y, report->_blockIdx.z);
      fprintf (stderr, "threadIdx=(%d %d %d)\n", report->_threadIdx.x, report->_threadIdx.y, report->_threadIdx.z);
    }
    if(file!=NULL)
      fprintf (stderr, "file %s, line %d\n", file, line);
    exit(1);
  }
}

#define CGBN_CHECK(report) cgbn_check(report, __FILE__, __LINE__)

// ---------------------------------------------------------------- //

// The CGBN context uses the following three parameters:
//   TPB             - threads per block (zero means to use the blockDim.x)
//   MAX_ROTATION    - must be small power of 2, imperically, 4 works well
//   CONSTANT_TIME   - require constant time algorithms (currently, constant time algorithms are not available)

// Locally it will also be helpful to have several parameters:
//   TPI             - threads per instance
//   BITS            - number of bits per instance

const uint32_t TPB_DEFAULT = 1024;

template<uint32_t tpi, uint32_t bits>
class cgbn_params_t {
  public:
  // parameters used by the CGBN context
  static const uint32_t TPB=TPB_DEFAULT;           // Reasonable default
  static const uint32_t MAX_ROTATION=4;            // good default value
  static const uint32_t SHM_LIMIT=0;               // no shared mem available
  static const bool     CONSTANT_TIME=false;       // not implemented

  // parameters used locally in the application
  static const uint32_t TPI=tpi;                   // threads per instance
  static const uint32_t BITS=bits;                 // instance size
};


// kernel implementation using cgbn
template<class params>
__global__ void kernel_iterate_cf(
        cgbn_error_report_t *report,
        uint64_t MAX_CF,
        uint32_t count,
        uint32_t *data,
        uint32_t *results
) {
  // decode an instance_i number from the blockIdx and threadIdx
  int32_t instance_i = (blockIdx.x*blockDim.x + threadIdx.x) / params::TPI;
  if(instance_i >= count)
    return;
  int32_t thread_i = (blockIdx.x*blockDim.x + threadIdx.x) % params::TPI;

  typedef cgbn_context_t<params::TPI, params>   context_t;
  typedef cgbn_env_t<context_t, params::BITS>   env_t;
  typedef typename env_t::cgbn_t                bn_t;
  typedef cgbn_mem_t<params::BITS>              mem_t;

  /* Cast uint32_t array to mem_t */
  mem_t *data_cast = (mem_t*) data;

  cgbn_monitor_t monitor = CHECK_ERROR ? cgbn_report_monitor : cgbn_no_checks;
  context_t _context(monitor, report, (uint32_t)instance_i);
  env_t _env(_context);

  bn_t x, a0, two_a0, a, b, c;
  bn_t t, t2;

  { // Setup
    cgbn_load(_env, x,  &data_cast[6*instance_i+0]);
    cgbn_load(_env, a0, &data_cast[6*instance_i+1]);
    cgbn_load(_env, two_a0, &data_cast[6*instance_i+2]);
    cgbn_load(_env, b,  &data_cast[6*instance_i+3]);
    cgbn_load(_env, c,  &data_cast[6*instance_i+4]);
    cgbn_load(_env, a,  &data_cast[6*instance_i+5]);
  }

  /*
  if (instance_i == 0 && thread_i == 0) {
    printf("GPU %2d | %u %u goal: %u\n", instance_i,
            cgbn_get_ui32(_env, x),
            cgbn_get_ui32(_env, a0),
            cgbn_get_ui32(_env, two_a0));
  }
  // */
  if (thread_i == 0)
    results[instance_i] = 0;

  for (uint64_t i = 2; i <= MAX_CF; i++) {
    // b = a*c - b;
    cgbn_mul(_env, t, a, c);
    cgbn_sub(_env, b, t, b);

    // c = (x - b*b) / c
    cgbn_sqr(_env, t, b);
    cgbn_sub(_env, t, x, t);

    // TODO: This didn't seem to improve timing much
    // test again, also consider testing for c = 2
    if (cgbn_equals_ui32(_env, c, 1)) {
        cgbn_set(_env, c, t); // Can skip this division because c == 1
        cgbn_add(_env, a, a0, b);
    } else {
        // 20% of runtime is this divide.
        cgbn_div(_env, c, t, c);
        // a = (a0 + b) / c
        cgbn_add(_env, t, a0, b);
        cgbn_div(_env, a, t, c); // DITTO
    }

    // Can this be moved under the c == 1 branch?
    // done = (a == two_a0)
    if (cgbn_equals(_env, two_a0, a) && results[instance_i] == 0) {
      results[instance_i] = i + 1;
    }
  }
}


void CgbnPessemisticCf::store_data(vector<pair<__uint128_t, __uint128_t>>& D_a0) {
    size_t    count = D_a0.size();
    assert( count <= max_D_count );
    for (size_t i = 0; i < count; i++) {
        auto [D, a0] = D_a0[i];
        auto c = D - a0 * a0;

        data[6*i + 0] = D;
        data[6*i + 1] = a0;
        data[6*i + 2] = a0 << 1; // two_a0
        data[6*i + 3] = a0; // b
        data[6*i + 4] = c; // c
        data[6*i + 5] = (a0 << 1) / c; // a
        //printf("%2lu | %lu %lu %lu\n", i, (uint64_t) D, (uint64_t) a0, (uint64_t) (a0 << 1));
    }
}


CgbnPessemisticCf::CgbnPessemisticCf(size_t D_size) {
    CUDA_CHECK(cudaEventCreate (&start));
    CUDA_CHECK(cudaEventCreate (&stop));

    max_D_count = D_size;
    data_size = 6 * sizeof(__uint128_t) * D_size;
    data = (__uint128_t*) malloc(data_size);

    CUDA_CHECK(cudaMalloc((void **)&gpu_data, data_size));
    CUDA_CHECK(cudaMalloc((void **)&gpu_results, sizeof(uint32_t) * D_size));
}

CgbnPessemisticCf::~CgbnPessemisticCf() {
    free(data);
    CUDA_CHECK(cudaFree(gpu_data));
    CUDA_CHECK(cudaFree(gpu_results));
    CUDA_CHECK(cudaEventDestroy (start));
    CUDA_CHECK(cudaEventDestroy (stop));
}

void CgbnPessemisticCf::run(size_t MAX_CF, vector<pair<__uint128_t, __uint128_t>>& D_a0, vector<uint32_t> &valid, int verbose) {
    // create a cgbn_error_report for CGBN to report back errors
    cgbn_error_report_t *report;
    CUDA_CHECK(cgbn_error_report_alloc(&report));

    size_t    count = D_a0.size();

    const int32_t   TPI =  4;
    const uint32_t  BITS = 128;      // kernel bits
    typedef cgbn_params_t<TPI, BITS>   cgbn_params_small;

    int32_t   TPB=TPB_DEFAULT; // Always the same default
    int32_t   IPB;             // IPB = TPB / TPI, instances per block
    size_t    BLOCK_COUNT;     // How many blocks to cover all count

    IPB = TPB / TPI;
    BLOCK_COUNT = (count + IPB - 1) / IPB;

    //kernel_info((const void*)kernel_iterate_cf<cgbn_params_small>, true);

    store_data(D_a0);

    CUDA_CHECK(cudaEventRecord (start));

    // Copy data
    if (verbose)
      printf("Copying %'lu bytes of data to GPU\n", data_size);
    CUDA_CHECK(cudaMemcpy(gpu_data, data, data_size, cudaMemcpyHostToDevice));

    if (verbose)
      printf("%lu -> CGBN<%d, %d> running kernel<%lu block x %d threads>\n",
          count, BITS, TPI, BLOCK_COUNT, TPB);

    /* Call CUDA Kernel. */
    kernel_iterate_cf<cgbn_params_small><<<BLOCK_COUNT, TPB>>>(
        report, MAX_CF, count, gpu_data, gpu_results);

    /* error report uses managed memory, sync the device and check for cgbn errors */
    CUDA_CHECK(cudaDeviceSynchronize());
    if (report->_error)
      printf("\n\nerror: %d\n", report->_error);
    CGBN_CHECK(report);

    float gputime;
    CUDA_CHECK(cudaEventRecord (stop));
    CUDA_CHECK(cudaEventSynchronize (stop));
    cudaEventElapsedTime (&gputime, start, stop);

    // Copy data back from GPU memory
    if (verbose)
      printf("Copying results back to CPU ...\n");
    CUDA_CHECK(cudaMemcpy(static_cast<void*>(valid.data()), gpu_results, sizeof(uint32_t) * count, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cgbn_error_report_free(report));
}

#ifdef __CUDA_ARCH__
  #if __CUDA_ARCH__ < 350
    #error "Unsupported architecture"
  #endif
#endif
