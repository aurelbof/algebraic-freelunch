#include <stdio.h>
#include <stdlib.h>
#include <flint/nmod_mpoly.h>

// #include "define.h"
#include "product.h"

// #if NVARS == 3
// #define DIMENSION_MATRIX 9
// #define DEGX 1029
// #elif NVARS == 4
// #define DIMENSION_MATRIX 27 // 3 ^ (R - 1) with R = 4
// #define DEGX 7203 // 3 * 5 * 7 ^ (R-3) with R = 4
// #elif NVARS == 5
// #define DIMENSION_MATRIX 81 // 3 ^ (R - 1) with R = 5
// #define DEGX 50421 // 3 * 5 * 7 ^ (R-3) with R = 5
// #elif NVARS == 6
// #define DIMENSION_MATRIX 243
// #define DEGX 36015
// #endif // NVARS == 3


// degree_array must have size nvars
void build_product(nmod_mpoly_struct *product, unsigned int *degree_array, nmod_mpoly_ctx_t ctx, const char **__x, long unsigned int *weights, unsigned int nvars, unsigned int degx)
{
    nmod_mpoly_struct **factors = (nmod_mpoly_struct **) malloc(nvars * sizeof(nmod_mpoly_struct *));
    nmod_mpoly_t den;
    nmod_mpoly_init(den, ctx);
    char strpol[50];
    memset(strpol, 0, 50);

    nmod_mpoly_one(product, ctx);
    for (unsigned int i = 0; i < nvars-1; ++i)
    {
        factors[i] = (nmod_mpoly_struct *) malloc(sizeof(nmod_mpoly_struct));
        nmod_mpoly_init(factors[i], ctx);
        snprintf(strpol, 50, "z%d^%ld-1", i+1, 5*weights[i]);
        nmod_mpoly_set_str_pretty(factors[i], strpol, __x, ctx);
        snprintf(strpol, 50, "z%d^%ld-1", i+1, weights[i]);
        nmod_mpoly_set_str_pretty(den, strpol, __x, ctx);
        nmod_mpoly_div(factors[i], factors[i], den, ctx);
        nmod_mpoly_print_pretty(factors[i], __x, ctx);
        puts("");
        nmod_mpoly_mul(product, product, factors[i], ctx);
    }

    // last factor
    factors[nvars-1] = (nmod_mpoly_struct *) malloc(sizeof(nmod_mpoly_struct));
    nmod_mpoly_init(factors[nvars-1], ctx);
    snprintf(strpol, 50, "z%d^%ld-1", nvars, (2*degx+1)*weights[nvars-1]);
    nmod_mpoly_set_str_pretty(factors[nvars-1], strpol, __x, ctx);
    snprintf(strpol, 50, "z%d^%ld-1", nvars, weights[nvars-1]);
    nmod_mpoly_set_str_pretty(den, strpol, __x, ctx);
    nmod_mpoly_div(factors[nvars-1], factors[nvars-1], den, ctx);
    nmod_mpoly_mul(product, product, factors[nvars-1], ctx);

    // nmod_mpoly_reverse(product, product, ctx);

    mp_limb_signed_t degrees[nvars];
    FILE *sorted_monomials = fopen("sorted_monomials.txt", "w");
    mp_limb_t N = mpoly_words_per_exp(product->bits, ctx->minfo);

    // uint degree_array[NVARS] = {5,5,5,1471};

    for (int i = 0; i < product->length; ++i)
    {
        mpoly_degrees_si(degrees, product->exps + N*i, 1, product->bits, ctx->minfo);
        //fprintf(sorted_monomials, "%ld ", multi_to_uni(degrees, degree_array));
    }
    fclose(sorted_monomials);

    nmod_mpoly_clear(den, ctx);
    for (unsigned int i = 0; i < nvars; ++i)
    {
        nmod_mpoly_clear(factors[i], ctx);
        free(factors[i]);
    }

    free(factors);
}

// int main()
// {
//     nmod_mpoly_ctx_t ctx;
//     nmod_mpoly_ctx_init(ctx, NVARS, ORD_DEGLEX, MODULUS);

//     nmod_mpoly_struct **factors = (nmod_mpoly_struct **) malloc(NVARS * sizeof(nmod_mpoly_struct *));
//     nmod_mpoly_t product, den;
//     nmod_mpoly_init(product, ctx);
//     nmod_mpoly_init(den, ctx);
//     char strpol[50];
//     memset(strpol, 0, 50);

//     nmod_mpoly_one(product, ctx);
//     for (int i = 0; i < NVARS-1; ++i)
//     {
//         factors[i] = (nmod_mpoly_struct *) malloc(sizeof(nmod_mpoly_struct));
//         nmod_mpoly_init(factors[i], ctx);
//         snprintf(strpol, 50, "z%d^%ld-1", i+1, 5*weights[i]);
//         nmod_mpoly_set_str_pretty(factors[i], strpol, __x, ctx);
//         snprintf(strpol, 50, "z%d^%ld-1", i+1, weights[i]);
//         nmod_mpoly_set_str_pretty(den, strpol, __x, ctx);
//         nmod_mpoly_div(factors[i], factors[i], den, ctx);
//         nmod_mpoly_print_pretty(factors[i], __x, ctx);
//         puts("");
//         nmod_mpoly_mul(product, product, factors[i], ctx);
//     }

//     // last factor
//     factors[NVARS-1] = (nmod_mpoly_struct *) malloc(sizeof(nmod_mpoly_struct));
//     nmod_mpoly_init(factors[NVARS-1], ctx);
//     snprintf(strpol, 50, "z%d^%ld-1", NVARS, (2*DEGX+1)*weights[NVARS-1]);
//     nmod_mpoly_set_str_pretty(factors[NVARS-1], strpol, __x, ctx);
//     snprintf(strpol, 50, "z%d^%ld-1", NVARS, weights[NVARS-1]);
//     nmod_mpoly_set_str_pretty(den, strpol, __x, ctx);
//     nmod_mpoly_div(factors[NVARS-1], factors[NVARS-1], den, ctx);
//     nmod_mpoly_mul(product, product, factors[NVARS-1], ctx);

//     // nmod_mpoly_reverse(product, product, ctx);

//     mp_limb_signed_t degrees[NVARS];
//     FILE *sorted_monomials = fopen("sorted_monomials.txt", "w");
//     mp_limb_t N = mpoly_words_per_exp(product->bits, ctx->minfo);

//     mp_limb_t degree_array[NVARS] = {5,5,5,1471};

//     for (int i = 0; i < product->length; ++i)
//     {
//         mpoly_degrees_si(degrees, product->exps + N*i, 1, product->bits, ctx->minfo);
//         fprintf(sorted_monomials, "%ld\n", multi_to_uni(degrees, degree_array));
//     }
//     fclose(sorted_monomials);

//     nmod_mpoly_clear(den, ctx);
//     nmod_mpoly_clear(product, ctx);
//     for (int i = 0; i < NVARS; ++i)
//     {
//         nmod_mpoly_clear(factors[i], ctx);
//         free(factors[i]);
//     }

//     free(factors); 
//     return 0;
// }
