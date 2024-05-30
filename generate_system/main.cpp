#include <iostream>
#include <NTL/lzz_pXFactoring.h>

#include <util.h>
#include <mat_lzz_pX_utils.h>
#include <mat_lzz_pX_determinant.h>

#include <flint/nmod_mpoly.h>
#include <flint/nmod_poly_mat.h>


#include "normal_form.h"
#include "generate_system.h"

#define COMPUTE_HARD_COLUMN 1

void decomp_multi_index(slong *indices, slong *degs, int c, int nvars);

int main(int argc, char **argv)
{
    if (argc != 2)
    {
        std::cerr << "Usage: ./main <nvars>" << std::endl;
        return EXIT_FAILURE;
    }
    int nvars = atoi(argv[1]);
    int tt = 12;
    long long MODULUS = 28407454060060787;
    nmod_mpoly_ctx_t ctx;
    nmod_mpoly_ctx_init(ctx, nvars, ORD_DEGLEX, MODULUS);

    NTL::zz_p::init(MODULUS);

    char** __x = (char **)malloc(nvars * sizeof(char*));
    printf("Hello");
    for (int i = 0; i < nvars; ++i)
    {
        __x[i] = (char*) malloc(8*sizeof(char));
        snprintf(__x[i], 8, "z%d", (ushort)(i+1));
    }
    printf("%s\n", __x[2]);
    static unsigned long int *weights = (unsigned long int *)malloc(nvars*sizeof(unsigned long int));
    int tmp = 7;
    for (int i = 0; i < nvars-1; ++i) {
        tmp *= 7;
        weights[i] = tmp+1;
    }
    weights[nvars-1] = 3;

    // first, generate the variables
    nmod_mpoly_t variables[nvars];
    char varstr[10];
    for (int i = 0; i < nvars; ++i)
    {
        nmod_mpoly_init(variables[i], ctx);
        snprintf(varstr, 10, "z%d^%ld", (short)(i+1), weights[i]);
        nmod_mpoly_set_str_pretty(variables[i], varstr, (const char **)__x, ctx);
    }

    puts("Beginning generation...");
    clock_t begin_generation = clock(), end_generation;
    nmod_mpoly_struct **F = generate_system(ctx, (const char **)__x, weights, nvars, tt);
    
    end_generation = clock();

    printf("System has been generated in %f s\n", ((double) (end_generation - begin_generation) / CLOCKS_PER_SEC));

    // Now write the data into a file
    char filename[100];
    snprintf(filename, 100, "../solve_system/equations_files/flint_system%d.txt", nvars);
    FILE *flint_system = fopen(filename, "w");
    fprintf(flint_system, "%d\n", nvars);
    for (int i = 0; i < nvars; ++i)
        fprintf(flint_system, "%s\n", __x[i]);

    for (int i = 0; i < nvars; ++i)
        fprintf(flint_system, "%ld\n", weights[i]);

    for (int i = 0; i < nvars; ++i)
    {
        nmod_mpoly_fprint_pretty(flint_system, F[i], (const char **)__x, ctx);
        fprintf(flint_system, "\n");
    }
    
    for (int i = 0; i < nvars; ++i)
    {
        nmod_mpoly_clear(F[i], ctx);
        free(F[i]);
    }
    free(F);

    free(__x);
    free(weights);
    return 0;
}