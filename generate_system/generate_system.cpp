#include <iostream>
#include <fstream>

#include <cstdlib>
#include <time.h>

// #include "define.h"
#include "generate_system.h"
#include "normal_form.h"
#include "fast_multivariate_multiplication.h"
#include "product.h"

#define STRLEN 1000

#define CLEAR(ff) nmod_mpoly_clear(ff, ctx)
#define FREE_SYSTEM(FF, ctx, nvars)          \
    for(int i = 0; i < nvars; ++i)    \
    {                                 \
        nmod_mpoly_clear(FF[i], ctx); \
        free(FF[i]);                  \
    }                                 \
    free(FF);       

using namespace NTL;

static mp_limb_t M3[3][3] = 
{
    {2, 1, 1},
    {1, 2, 1},
    {1, 1, 2}
};

static mp_limb_t M4[4][4] = 
{
    {5, 7, 1, 3},
    {4, 6, 1, 1},
    {1, 3, 5, 7},
    {1, 1, 4, 6}
};

/**
 * @brief build the matrix used in the linear layer of Griffin
*/
mp_limb_t ** build_matrix(int tt)
{   
    mp_limb_t **MDS_LAYER = (mp_limb_t **)malloc(tt * sizeof(mp_limb_t *));
    for (int i = 0; i < tt; ++i)
        MDS_LAYER[i] = (mp_limb_t *)malloc(tt * sizeof(mp_limb_t));
    if (tt == 3)
    {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                MDS_LAYER[i][j] = M3[i][j];
    }

    if (tt == 4)
    {
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                MDS_LAYER[i][j] = M4[i][j];
    }

    else
    {
        int tp = tt/4;
        for (int k = 0; k < tp; ++k)
        {
            for (int l = 0; l < tp; ++l)
            {
                mp_limb_t c = k == l ? 2 : 1;
                for (int i = 0; i < 4; ++i)
                    for (int j = 0; j < 4; ++j)
                        MDS_LAYER[4*k+i][4*l+j] = c * M4[i][j]; 
            }
        }
    }
    return MDS_LAYER;
}


// void fill_first_round(nmod_mpoly_struct *P, nmod_mpoly_struct **Q0, nmod_mpoly_t *variables, nmod_mpoly_ctx_t ctx, const char **__x, long unsigned int *weights, int nvars, int tt);

nmod_mpoly_struct **generate_system(nmod_mpoly_ctx_struct *ctx, const char **__x, long unsigned int *weights, int nvars, int tt)
{
    int degx = 3 * 5;
    for (int _ = 0; _ < nvars - 1; ++_)
        degx *= 7;
    printf("degx = %d\n", degx);
    
    long long MODULUS = ctx->mod.n;

    mp_limb_t **MDS_LAYER = build_matrix(tt);
    
    for (int i = 0; i < tt; ++i)
    {
        for (int j = 0; j < tt; ++j)
            printf("%ld ", MDS_LAYER[i][j]);
        printf("\n");
    }

    // build the polynomial "product"
    nmod_mpoly_t product;
    nmod_mpoly_init(product, ctx);
    uint degree_array[nvars];
    for (int i = 0; i < nvars-1; ++i)
        degree_array[i] = 5;
    degree_array[nvars-1] = 21;
    for (int i = 0; i < nvars-1; ++i)
        degree_array[nvars-1] *= 7;
    degree_array[nvars-1] = 2*degree_array[nvars-1]+1;

    // debug
    for (int i = 0; i < nvars; ++i)
        printf("degree_array[%d] = %d\n", i, degree_array[i]);
    puts("");

    // test degree_array
    for (int i = 0; i < nvars; ++i)
        printf("%d   ", degree_array[i]);

    // zz_p::init(MODULUS);

    puts("building the product...");
    build_product(product, degree_array, ctx, __x, weights, nvars, degx);
    puts("Done.");

    // Read the round constants in the appropriate file
    mp_limb_t constants[tt * nvars];
    char file_constants_name[50];
    snprintf(file_constants_name, 50, "constants%d.txt", nvars);
    std::ifstream file_constants(file_constants_name);

    for (int i = 0; i < tt*nvars; ++i)
        file_constants >> constants[i];

    
    FILE *txtproduct = fopen("product.txt", "w");
    nmod_mpoly_fprint_pretty(txtproduct, product, __x, ctx);
    fclose(txtproduct);

    nmod_mpoly_t Psquare, tmp1, tmp2;
    nmod_mpoly_init(Psquare, ctx);
    nmod_mpoly_init(tmp1, ctx);
    nmod_mpoly_init(tmp2, ctx);
    nmod_mpoly_struct **F = (nmod_mpoly_struct **) malloc(nvars * sizeof(nmod_mpoly_struct *));
    nmod_mpoly_struct **P = (nmod_mpoly_struct **) malloc(nvars * sizeof(nmod_mpoly_struct *));
    nmod_mpoly_struct ***Q = (nmod_mpoly_struct ***) malloc(nvars * sizeof(nmod_mpoly_struct **));
    
    for (int i = 0; i < nvars; ++i)
    {
        F[i] = (nmod_mpoly_struct *) malloc(sizeof(nmod_mpoly_struct));
        P[i] = (nmod_mpoly_struct *) malloc(sizeof(nmod_mpoly_struct));

        nmod_mpoly_init(F[i], ctx);
        nmod_mpoly_init(P[i], ctx);
    }

    // initialization of Q
    for (int i = 0; i < nvars; ++i)
    {
        Q[i] = (nmod_mpoly_struct **) malloc((tt-2) * sizeof(nmod_mpoly_struct*));
        for (int j = 0; j < tt-2; ++j)
        {
            Q[i][j] = (nmod_mpoly_struct *) malloc(sizeof(nmod_mpoly_struct));
            nmod_mpoly_init(Q[i][j], ctx);
        }
    }

    // build an array containing each variable
    nmod_mpoly_t variables[nvars+1];
    char varstr[10];
    nmod_mpoly_init(variables[0], ctx);
    nmod_mpoly_one(variables[0], ctx);
    nmod_t N;
    nmod_init(&N, MODULUS);
    for (int i = 1; i <= nvars; ++i)
    {
        nmod_mpoly_init(variables[i], ctx);
        snprintf(varstr, 10, "z%d^%ld", (short)i, weights[i-1]);
        nmod_mpoly_set_str_pretty(variables[i], varstr, __x, ctx);
    }
    for (int i = 0; i <= nvars; ++i)
    {
        printf("DEBUG: variables[%d] = ", i);
        nmod_mpoly_print_pretty(variables[i], __x, ctx);
        puts("");
    }

    /* Before the first round */
    nmod_mpoly_t sum_right, lhs, rhs;
    nmod_mpoly_init(sum_right, ctx);
    nmod_mpoly_init(rhs, ctx);
    nmod_mpoly_init(lhs, ctx);

    mp_limb_t alpha_i = 4, beta_i = 7, gamma_i = 1;

    char filename[40];
    memset(filename, 0, 40);
    snprintf(filename, 40, "third_state.txt");
    FILE *file = fopen(filename, "r");

    char *strpoly = (char *) malloc(STRLEN * sizeof(char));
    memset(strpoly, 0, STRLEN);
    if (fscanf(file, "%s\n", strpoly) != 1)
    {
        fprintf(stderr, "Error while reading %s (1)", filename);
        exit(EXIT_FAILURE);
    }
    nmod_mpoly_set_str_pretty(variables[0], strpoly, __x, ctx);

    memset(strpoly, 0, STRLEN);
    if (fscanf(file, "%s\n", strpoly) != 1)
    {
        fprintf(stderr, "Error while reading %s (2)", filename);
        exit(EXIT_FAILURE);
    }
    nmod_mpoly_set_str_pretty(P[0], strpoly, __x, ctx);

    for (int j = 0; j < tt-2; ++j)
    {
        memset(strpoly, 0, STRLEN);
        if (fscanf(file, "%s\n", strpoly) != 1)
        {
            fprintf(stderr, "Error while reading %s (%d)", filename, j+3);
            exit(EXIT_FAILURE);
        }
        nmod_mpoly_set_str_pretty(Q[0][j], strpoly, __x, ctx);
    }

    fclose(file);
    free(strpoly);

    printf("var[0] = ");
    nmod_mpoly_print_pretty(variables[0], __x, ctx);
    puts("");

    printf("P0 = ");
    nmod_mpoly_print_pretty(P[0], __x, ctx);
    puts("");

    for (int j = 0; j < tt-2; ++j)
    {
        printf("Q0[%d] = ", j);
        nmod_mpoly_print_pretty(Q[0][j], __x, ctx);
        puts("");
    }

    clock_t equation_start, equation_end, qi_start, qi_end;
    double cpu_time_used;

    equation_start = clock();
    // just for initialization
    equation_end = clock();
    printf("Generating equation %d\n", 0);
    fflush(stdout);
    for (int i = 0; i < nvars-1; ++i)
    {
        nmod_mpoly_pow_ui(F[i], variables[i+1], 3, ctx);
        nmod_mpoly_scalar_addmul_ui(F[i], F[i], variables[i], MODULUS - MDS_LAYER[0][0], ctx);
        nmod_mpoly_scalar_addmul_ui(F[i], F[i], P[i], MODULUS - MDS_LAYER[0][1], ctx);
        for (int j = 0; j < tt-2; ++j)
            nmod_mpoly_scalar_addmul_ui(F[i], F[i], Q[i][j], MODULUS - MDS_LAYER[0][j+2], ctx);

        nmod_mpoly_sub_ui(F[i], F[i], constants[tt*i], ctx);  

        equation_end = clock();
        cpu_time_used = ((double) (equation_end - equation_start)) / CLOCKS_PER_SEC;
        printf("TIME FOR EQUATION %d: %f s\n", i, cpu_time_used);
        fflush(stdout);

        equation_start = clock();
        printf("Generating equation %d\n", i+1);
        fflush(stdout);
        // nmod_mpoly_scalar_addmul_ui(A, B, C, c, ctx): set A to B + c*C

        // first compute P[i+1] = (m_{1,0}, m_{1,1}, m_{1,2}, ... , m_{1,tt-1}) . (zi, P[i], Q[i][0], Q[i][1], ..., Q[i][tt-3]) + cste

        nmod_mpoly_scalar_mul_ui(P[i+1], variables[i], MDS_LAYER[1][0], ctx);
        nmod_mpoly_scalar_addmul_ui(P[i+1], P[i+1], P[i], MDS_LAYER[1][1], ctx);

        for (int j = 0; j < tt-2; ++j)
            nmod_mpoly_scalar_addmul_ui(P[i+1], P[i+1], Q[i][j], MDS_LAYER[1][j+2], ctx);

        nmod_mpoly_add_ui(P[i+1], P[i+1], constants[tt*i + 1], ctx);

        puts("Computing P^3");
        fflush(stdout);
        clock_t pcube = clock();
        // nmod_mpoly_pow_ui(P[i+1], P[i+1], 3, ctx);
        if (i < 2)
            nmod_mpoly_pow_ui(Psquare, P[i+1], 2, ctx);
        else
            multiply(P[i+1], P[i+1], Psquare, degree_array, product, ctx, weights, nvars);
        normal_form(Psquare, F, i , ctx);
        if (i < 2)
            nmod_mpoly_mul(P[i+1], Psquare, P[i+1], ctx);
        else
            multiply(P[i+1], Psquare, P[i+1], degree_array, product, ctx, weights, nvars);
        normal_form(P[i+1], F, i, ctx);

        printf("P^3 computed in %f s\n", (double) (clock() - pcube) / CLOCKS_PER_SEC);

        fflush(stdout);

        // compute P^2 and reduce
        if (i < 2)
            nmod_mpoly_mul(Psquare, P[i+1], P[i+1], ctx);
        else
            multiply(P[i+1], P[i+1], Psquare, degree_array, product, ctx, weights, nvars);
        normal_form(Psquare, F, i, ctx);

        // now compute each polynomial Q[i+1][j]
        alpha_i = 4; beta_i = 7; gamma_i = 1;
        for (int j = 0; j < tt-2; ++j)
        {
            qi_start = clock();
            // as above, first compute the rhs
            nmod_mpoly_scalar_addmul_ui(sum_right, P[i+1], variables[i+1], gamma_i, ctx);
            if (j > 0)
                nmod_mpoly_add(sum_right, sum_right, lhs, ctx);
            else
                nmod_mpoly_zero(lhs, ctx);

            // nmod_mpoly_mul(rhs, sum_right, sum_right, ctx);                   // rhs = P[0]^2
            // nmod_mpoly_scalar_addmul_ui(rhs, rhs, sum_right, alpha_i, ctx);     // rhs = P[0]^2 + a_i*P[0]
            // nmod_mpoly_add_ui(rhs, rhs, beta_i, ctx);
            
            nmod_mpoly_scalar_addmul_ui(tmp1, lhs, variables[i+1], gamma_i, ctx);

            // compute tmp2 right away
            if (i < 2)
                nmod_mpoly_pow_ui(tmp2, tmp1, 2, ctx);
            else
                multiply(tmp1,tmp1,tmp2,degree_array,product,ctx,weights,nvars);
            // Augustin: I changed the previous lines: a bug might come from here.
            normal_form(tmp2, F, i, ctx);
            if (i < 2)
                nmod_mpoly_mul(tmp1, tmp1, P[i+1], ctx);
            else
                multiply(tmp1, P[i+1], tmp1, degree_array, product, ctx, weights, nvars);
            
            normal_form(tmp1, F, i, ctx);
            nmod_mpoly_scalar_mul_ui(tmp1, tmp1, 2, ctx); // tmp1 = 2*P[i+1] * (gamma Z_{i+1} + lhs') mod I

            nmod_mpoly_add(rhs, tmp1, tmp2, ctx);
            nmod_mpoly_add(rhs, rhs, Psquare, ctx);

            nmod_mpoly_scalar_addmul_ui(rhs, rhs, sum_right, alpha_i, ctx);
            nmod_mpoly_add_ui(rhs, rhs, beta_i, ctx);

            // now compute the current lhs
            nmod_mpoly_scalar_mul_ui(lhs, variables[i], MDS_LAYER[j+2][0], ctx);
            nmod_mpoly_scalar_addmul_ui(lhs, lhs, P[i], MDS_LAYER[j+2][1], ctx);

            for (int k = 0; k < tt-2; ++k)
                nmod_mpoly_scalar_addmul_ui(lhs, lhs, Q[i][k], MDS_LAYER[j+2][k+2], ctx);

            nmod_mpoly_add_ui(lhs, lhs, constants[tt*i+j+2], ctx);
            if (i < 2)
                nmod_mpoly_mul(Q[i+1][j], lhs, rhs, ctx);
            else
                multiply(lhs, rhs, Q[i+1][j], degree_array, product, ctx, weights, nvars);                

            clock_t nf_Q_start = clock();
            normal_form(Q[i+1][j], F, i, ctx);
            printf("Normal form of Q computed in %f s\n", (double) (clock() - nf_Q_start) / CLOCKS_PER_SEC);

            fflush(stdout);
            qi_end = clock();
            if (i >= 2)
                printf("Q[%d][%d] has been generated in %f s\n", i, j, ((double)(qi_end - qi_start) / CLOCKS_PER_SEC));

            fflush(stdout);
            alpha_i = (j+2) * 4;
            beta_i  = (j+2)*(j+2) * 7;
            gamma_i++;
        }        
    }

    // don't forget the last equation
    
    printf("Computing the last equation\n");
    fflush(stdout);

    nmod_mpoly_scalar_mul_ui(F[nvars-1], variables[nvars-1], MDS_LAYER[0][0], ctx);
    nmod_mpoly_scalar_addmul_ui(F[nvars-1], F[nvars-1], P[nvars-1], MDS_LAYER[0][1], ctx);
    for (int j = 0; j < tt-2; ++j)
        nmod_mpoly_scalar_addmul_ui(F[nvars-1], F[nvars-1], Q[nvars-1][j], MDS_LAYER[0][j+2], ctx);

    nmod_mpoly_make_monic(F[nvars-1], F[nvars-1], ctx);

    equation_start = clock();
    cpu_time_used = ((double) (equation_start - equation_end)) / CLOCKS_PER_SEC;
    printf("TIME FOR LAST EQUATION: %f s\n", cpu_time_used);
    fflush(stdout);
    
    for (int i = 0; i < tt; ++i) {
        free(MDS_LAYER[i]);
    }
    free(MDS_LAYER);

    CLEAR(product);
    CLEAR(Psquare);
    CLEAR(tmp1);
    CLEAR(tmp2);
    CLEAR(lhs);
    CLEAR(rhs);
    CLEAR(sum_right);
    FREE_SYSTEM(P, ctx, nvars);


    
    // free Q
    for (int i = 0; i < nvars; ++i)
    {
        for (int j = 0; j < tt-2; ++j)
        {
            nmod_mpoly_clear(Q[i][j], ctx);
            free(Q[i][j]);
        }
        free(Q[i]);
    }
    free(Q);

    for (int i = 0; i <= nvars; ++i)
        nmod_mpoly_clear(variables[i], ctx);
    return F;
}