#ifndef PRODUCT_H
#define PRODUCT_H

#include <sys/types.h>
#include <flint/nmod_mpoly.h>

static inline unsigned long multi_to_uni(mp_limb_signed_t *degs, unsigned int *degree_array, long unsigned int *weights, int nvars)
{
    unsigned long curr_var_deg = 1;
    unsigned long monom_eq_deg = 0;
    for (int j = 0; j < nvars; ++j)
    {
        monom_eq_deg += curr_var_deg*degs[j]/weights[j];
        curr_var_deg *= (ulong) degree_array[j];
    }
    return monom_eq_deg;
}

void build_product(nmod_mpoly_struct *product, unsigned int *degree_array, nmod_mpoly_ctx_t ctx, const char **__x, long unsigned int *weights, unsigned int nvars, unsigned int degx);

#endif // PRODUCT_H
