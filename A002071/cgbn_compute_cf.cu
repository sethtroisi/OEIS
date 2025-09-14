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
//   TBP             - threads per block (zero means to use the blockDim.x)
//   MAX_ROTATION    - must be small power of 2, imperically, 4 works well
//   CONSTANT_TIME   - require constant time algorithms (currently, constant time algorithms are not available)

// Locally it will also be helpful to have several parameters:
//   TPI             - threads per instance
//   BITS            - number of bits per instance

/* TODO test how this changes gpu_throughput_test */
/* NOTE: >= 512 may not be supported for > 2048 bit kernels */
const uint32_t TPB_DEFAULT = 256;

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
        uint64_t i_start,
        uint64_t i_end,
        uint32_t *data,
        uint32_t count
        ) {
  // decode an instance_i number from the blockIdx and threadIdx
  int32_t instance_i = (blockIdx.x*blockDim.x + threadIdx.x)/params::TPI;
  if(instance_i >= count)
    return;

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
  bn_t t, equal;

  { // Setup
    cgbn_load(_env, x,  &data_cast[6*instance_i+0]);
    cgbn_load(_env, a0, &data_cast[6*instance_i+1]);
    cgbn_load(_env, two_a0, &data_cast[6*instance_i+2]);
    cgbn_load(_env, a,  &data_cast[6*instance_i+3]);
    cgbn_load(_env, b,  &data_cast[6*instance_i+4]);
    cgbn_load(_env, c,  &data_cast[6*instance_i+5]);
  }

  data[instance_i] = 0;

  for (uint64_t i = 2; i < i_end; i++) {
    // b = a*c - b;
    cgbn_mul(_env, t, a, c);
    cgbn_sub(_env, b, t, b);

    // c = (x - b*b) / c
    cgbn_sqr(_env, t, b);
    cgbn_sub(_env, t, x, t);
    cgbn_div(_env, c, t, c);

    // a = (a0 + b) / c
    cgbn_add(_env, t, a0, b);
    cgbn_div(_env, a, t, c);

    // done = (a == two_a0)
    if (cgbn_equals(_env, a, two_a0) == 0) {
      if (data[instance_i] == 0) {
        data[instance_i] = i;
      }
    }
  }
}



void cgbn_pessemistic_cf(vector<__uint128_t>& D, vector<uint32_t> &valid, int verbose)
{

  cudaEvent_t start, stop;
  CUDA_CHECK(cudaEventCreate (&start));
  CUDA_CHECK(cudaEventCreate (&stop));
  CUDA_CHECK(cudaEventRecord (start));

  cgbn_error_report_t *report;
  // create a cgbn_error_report for CGBN to report back errors
  CUDA_CHECK(cgbn_error_report_alloc(&report));

  size_t    count = D.size();
  size_t    data_size;
  uint32_t  *gpu_data;

  uint32_t  BITS = 0;        // kernel bits
  int32_t   TPB=TPB_DEFAULT; // Always the same default
  int32_t   TPI;
  int32_t   IPB;             // IPB = TPB / TPI, instances per block
  size_t    BLOCK_COUNT;     // How many blocks to cover all count

  typedef cgbn_params_t<4, 256>   cgbn_params_small;
  TPI = cgbn_params_small::TPI;
  IPB = TPB / TPI;
  BLOCK_COUNT = (count + IPB - 1) / IPB;

  //kernel_info((const void*)kernel_iterate_cf<cgbn_params_small>, verbose);

  /* Consistency check that struct cgbn_mem_t is byte aligned without extra fields. */
  data_size = sizeof(__uint128_t) * D.size();

  // Copy data
  printf("Copying %'lu bytes of data to GPU\n", data_size);
  CUDA_CHECK(cudaMalloc((void **)&gpu_data, data_size));
  CUDA_CHECK(cudaMemcpy(gpu_data, static_cast<void*>(D.data()), data_size, cudaMemcpyHostToDevice));

  printf("CGBN<%d, %d> running kernel<%lu block x %d threads>\n",
      BITS, TPI, BLOCK_COUNT, TPB);

  /* Call CUDA Kernel. */
  kernel_iterate_cf<cgbn_params_small><<<BLOCK_COUNT, TPB>>>(
      report, 0, 200, gpu_data, count);

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
  printf("Copying results back to CPU ...\n");
  CUDA_CHECK(cudaMemcpy(static_cast<void*>(valid.data()), gpu_data, data_size, cudaMemcpyDeviceToHost));

  // clean up
  CUDA_CHECK(cudaFree(gpu_data));
  CUDA_CHECK(cgbn_error_report_free(report));
  CUDA_CHECK(cudaEventDestroy (start));
  CUDA_CHECK(cudaEventDestroy (stop));
  return;
}

#ifdef __CUDA_ARCH__
  #if __CUDA_ARCH__ < 350
    #error "Unsupported architecture"
  #endif
#endif
