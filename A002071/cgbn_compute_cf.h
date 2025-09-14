/* cgbn_compute_cf.h: header for CGBN (GPU) based continued fractions.

  Copyright 2025 Seth Troisi

  This program is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 3 of the License, or (at your
  option) any later version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
  more details.

  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
  02111-1307, USA.
*/

#pragma once

#include <cstdint>
#include <vector>
#include <utility>

// GMP import must proceed cgbn.h
#include <cuda_runtime.h>

using std::pair;
using std::vector;

extern "C" {

class CgbnPessemisticCf
{
    private:
        cudaEvent_t start, stop;
        size_t max_D_count;
        size_t    data_size;
        __uint128_t *data;
        uint32_t  *gpu_data;
        uint32_t  *gpu_results;

        void store_data(vector<pair<__uint128_t, __uint128_t>>& D_a0);

    public:
        CgbnPessemisticCf(size_t D_size);
        ~CgbnPessemisticCf();
        CgbnPessemisticCf(const CgbnPessemisticCf& that) = delete;
        CgbnPessemisticCf& operator=(const CgbnPessemisticCf& that) = delete;
        void run(size_t MAX_CF, vector<pair<__uint128_t, __uint128_t>>& D_a0, vector<uint32_t> &valid, int verbose);
};

}
