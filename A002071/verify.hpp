#pragma once

#include <vector>
#include <gmpxx.h>

using std::pair;
using std::vector;

vector<mpz_class> power_set_verify(const vector<uint32_t>& set);
pair<mpz_class, mpz_class> pell_PQA_verify(mpz_class D);

// Do a single iteration of pell_solution_cf(D)
void verify_expand_D(char* argv1);
