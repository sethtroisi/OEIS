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

using std::pair;
using std::vector;

extern "C" {

void cgbn_pessemistic_cf(size_t MAX_CF, vector<pair<__uint128_t, __uint128_t>>& D, vector<uint32_t> &valid, int verbose);

}
