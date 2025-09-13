#pragma once

#include <cstdint>
#include <vector>
#include <gmpxx.h>

using std::pair;
using std::vector;

vector<mpz_class> power_set_verify(const vector<uint32_t>& set);

int test_smooth_small_verify(mpz_class n, vector<uint32_t>& primes);

pair<mpz_class, mpz_class> pell_PQA_verify(mpz_class D);
