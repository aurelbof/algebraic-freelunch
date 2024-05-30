#ifndef FAST_MULTIVARIATE_MULTIPLICATION_H
#define FAST_MULTIVARIATE_MULTIPLICATION_H

#include <flint/nmod_mpoly.h>
#include <NTL/lzz_pXFactoring.h>
#include <NTL/ZZ_pX.h>

// #include "define.h"

void ntl_uni_to_flint_multi(nmod_mpoly_t P, NTL::zz_pX &old_pol, uint *degree_array, nmod_mpoly_struct const * product, nmod_mpoly_ctx_struct *ctx, long unsigned int *weights, int nvars);

void flint_multi_to_ntl_uni(nmod_mpoly_t const P, NTL::zz_pX &new_pol, uint *degree_array, nmod_mpoly_ctx_struct *ctx,long unsigned int *weights, int nvars);

void multiply(nmod_mpoly_t const P, nmod_mpoly_t const Q, nmod_mpoly_t PQ, uint *degree_array, nmod_mpoly_struct const * product, nmod_mpoly_ctx_struct *ctx, long unsigned int *weights, int nvars);

#endif // FAST_MULTIVARIATE_MULTIPLICATION_H
