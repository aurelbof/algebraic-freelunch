#ifndef GENERATE_SYSTEM_H
#define GENERATE_SYSTEM_H

#include <flint/nmod_mpoly.h>

// #include "define.h"

/**
 * @brief build the algrebraic system
 * @return an array of polynomials (f1, ..., f_{R-1}, g)
*/
nmod_mpoly_struct **generate_system(nmod_mpoly_ctx_struct *ctx, const char **__x, long unsigned int *weights, int nvars, int tt);

#endif // GENERATE_SYSTEM_H
