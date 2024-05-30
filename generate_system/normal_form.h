#ifndef NORMAL_FORM_H
#define NORMAL_FORM_H

#include <flint/nmod_mpoly.h>

/**
 * @brief set h to h mod F where F = [F[0], ... , F[N-1]]
 */
void normal_form(nmod_mpoly_t h, nmod_mpoly_struct **F, int N, nmod_mpoly_ctx_t ctx);

#endif // NORMAL_FORM_H
