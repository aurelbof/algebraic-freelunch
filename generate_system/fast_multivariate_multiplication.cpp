#include <stdio.h>
#include <time.h>

#include "fast_multivariate_multiplication.h"
// #include "read_poly.h"
#include "product.h"

using namespace NTL;

// Not implemented yet
void ntl_uni_to_flint_multi(nmod_mpoly_t P, zz_pX &old_pol, uint *degree_array, nmod_mpoly_struct const * product, nmod_mpoly_ctx_struct *ctx, long unsigned int *weights, int nvars)
{
    // printf("ntl -> flint\n");
    nmod_mpoly_zero(P, ctx);

    // clock_t begin = clock();
    mp_limb_signed_t degrees[nvars];
    uint _deg;
    mp_limb_t N = mpoly_words_per_exp(product->bits, ctx->minfo);
    for (int i = 0; i < product->length; ++i)
    {
        mpoly_degrees_si(degrees, product->exps + N*i, 1, product->bits, ctx->minfo);
        _deg = multi_to_uni(degrees, degree_array, weights, nvars);
        if (_deg <= deg(old_pol))
        {
            mp_limb_t coeff = conv<mp_limb_t, zz_p>(old_pol[_deg]);
            if (coeff != 0)
                nmod_mpoly_set_coeff_ui_ui(P, coeff, (const mp_limb_t *) degrees, ctx);
        }
    }
    // printf("Time ntl -> flint: %f s\n", (double) (clock() - begin) / CLOCKS_PER_SEC);
}


void flint_multi_to_ntl_uni(nmod_mpoly_t const P, zz_pX &new_pol, uint *degree_array, nmod_mpoly_ctx_struct *ctx,long unsigned int *weights, int nvars)
{
    // clock_t beg = clock();
    // printf("flint -> ntl\n");
    mp_limb_signed_t degs[nvars];
    mp_limb_t N = mpoly_words_per_exp(P->bits, ctx->minfo);
    for (int i = 0; i < P->length; ++i)
    {
        mpoly_degrees_si(degs, P->exps + N*i, 1, P->bits, ctx->minfo);

        unsigned long curr_var_deg = 1;
        unsigned long monom_eq_deg = 0;
        zz_p coeff(P->coeffs[i]);

        for (int j = 0; j < nvars; ++j)
        {
            monom_eq_deg += curr_var_deg*degs[j]/weights[j];
            curr_var_deg *= (ulong) degree_array[j];
        }
        SetCoeff(new_pol, monom_eq_deg, coeff);
    }
    // printf("Time for flint->ntl: %f\n", (double) (clock() - beg) / CLOCKS_PER_SEC);
}

// {5,5,5,5,5,7^{2R}}
// Typically degree_array = {5,5,1,2,2*d}
// z0 => X
// z1 => z0^5 = X^5
// z2 => z1^5 = X^{25}
// z3 => z2^1 
// z0^2 * z1^3 * z4^10 => (X)^2 * (X^5)^3 * (X^25)^{2*d}
// m = z0^a z1^b ... => X^{d_uni}

void multiply(nmod_mpoly_t const P, nmod_mpoly_t const Q, nmod_mpoly_t PQ, uint *degree_array, nmod_mpoly_struct const * product, nmod_mpoly_ctx_struct *ctx, long unsigned int *weights, int nvars)
{
    zz_pX P_NTL;
    flint_multi_to_ntl_uni(P, P_NTL, degree_array, ctx, weights, nvars);
    zz_pX Q_NTL;
    flint_multi_to_ntl_uni(Q, Q_NTL, degree_array,ctx, weights, nvars);

    zz_pX PQ_NTL;
    // clock_t beg = clock();
    // printf("Trying to multiply NTL polynomials of degree %i and %i \n", deg(P_NTL),deg(Q_NTL));
    // fflush(stdout);
    mul(PQ_NTL, P_NTL, Q_NTL);
    // printf("Time for mul: %f\n", (double) (clock() - beg) / CLOCKS_PER_SEC);

    ntl_uni_to_flint_multi(PQ, PQ_NTL, degree_array, product, ctx, weights, nvars);
}
