#include <iostream>
#include <NTL/lzz_pXFactoring.h>

#include <util.h>
#include <mat_lzz_pX_utils.h>
#include <mat_lzz_pX_determinant.h>

#include <flint/nmod_mpoly.h>
#include <flint/nmod_poly_mat.h>

#include "normal_form.h"

#define STRLEN 10000000

#define COMPUTE_HARD_COLUMN 1

void decomp_multi_index(slong *indices, slong *degs, int c, int nvars);

int main(int argc, char **argv)
{
    if (argc != 2)
    {
        fprintf(stderr, "Usage, ./main <path to equations file>\n");
        return EXIT_FAILURE;
    }
    int tt = 12;
    long long MODULUS = 28407454060060787;

    FILE *file = fopen(argv[1], "r");
    if (file == NULL)
    {
        fprintf(stderr, "Error: cannot open requested file.\n");
        exit(EXIT_FAILURE);
    }

    int nvars;
    fscanf(file, "%d\n", &nvars);

    nmod_mpoly_ctx_t ctx;
    nmod_mpoly_ctx_init(ctx, nvars, ORD_DEGLEX, MODULUS);

    clock_t begin_generation = clock(), end_generation;

    char **__x = (char**) malloc(nvars * sizeof(char*));

    for (int i = 0; i < nvars; ++i)
    {
        __x[i] = (char*)malloc(5*sizeof(char));
        fscanf(file, "%s\n", __x[i]);
    }


    std::cout << "nvars = " << nvars << std::endl;

    int weights[nvars];
    for(int i = 0; i < nvars; ++i)
        fscanf(file, "%d\n", weights+i);

    // now generate the variables
    nmod_mpoly_t variables[nvars];
    char varstr[10];
    for (int i = 0; i < nvars; ++i)
    {
        nmod_mpoly_init(variables[i], ctx);
        snprintf(varstr, 10, "%s^%d", __x[i], weights[i]);
        printf("varstr = %s\n", varstr);
        nmod_mpoly_set_str_pretty(variables[i], varstr, (const char**)__x, ctx);
    }

    for (int i = 0; i < nvars; ++i)
    {
        nmod_mpoly_print_pretty(variables[i], (const char **)__x, ctx);
        puts("");
    }
    nmod_mpoly_struct **F = (nmod_mpoly_struct **) malloc(nvars*sizeof(nmod_mpoly_struct*));

    for (int i = 0; i < nvars; ++i)
    {
        F[i] = (nmod_mpoly_struct *) malloc(sizeof(nmod_mpoly_struct));
        nmod_mpoly_init(F[i], ctx);
    }

    char *strpoly = (char *) malloc(STRLEN * sizeof(char));

    for (int i = 0; i < nvars; ++i)
    {
        memset(strpoly, 0, STRLEN);
        fscanf(file, "%s\n", strpoly);
        nmod_mpoly_set_str_pretty(F[i], strpoly, (const char **)__x, ctx);
    }

    free(strpoly);
    fclose(file);

    end_generation = clock();

    printf("System has been generated in %f s\n", ((double) (end_generation - begin_generation) / CLOCKS_PER_SEC));

    /* COMPUTE THE DIMENSION OF THE MATRIX */
    int dimension_matrix = 1;
    slong tmp_degs[nvars], degs[nvars];
    for (int i = 0; i < nvars-1; ++i)
    {
        // compute the degree of F[i] in zi
        nmod_mpoly_degrees_si(tmp_degs, F[i], ctx);
        degs[i] = tmp_degs[i] / weights[i];
        dimension_matrix *= degs[i];
    }

    nmod_mpoly_degrees_si(tmp_degs, F[nvars-1], ctx);
    degs[nvars-1] = tmp_degs[nvars-1] / weights[nvars-1];
    slong degx = degs[nvars-1];

    nmod_mpoly_t *monomials_to_reduce = (nmod_mpoly_t *) malloc(dimension_matrix * sizeof(nmod_mpoly_t));
    nmod_mpoly_t *monomial_basis = (nmod_mpoly_t *) malloc(dimension_matrix * sizeof(nmod_mpoly_t));

    char str_monomial[100];

    nmod_mpoly_t xdegx;
    nmod_mpoly_init(xdegx, ctx);
    nmod_mpoly_pow_ui(xdegx, variables[nvars-1], degx, ctx);


    // create the monomial basis and the monomials to be reduced
    slong exponents[nvars];
    nmod_mpoly_t _var;
    nmod_mpoly_init(_var, ctx);
    for (int c = 0; c < dimension_matrix; ++c)
    {
        nmod_mpoly_init(monomials_to_reduce[c], ctx);
        nmod_mpoly_init(monomial_basis[c], ctx);
        decomp_multi_index(exponents, degs, c, nvars);
        nmod_mpoly_one(monomial_basis[c], ctx);
        
        for (int i = 0; i < nvars; ++i)
        {
            nmod_mpoly_pow_ui(_var, variables[i], exponents[i], ctx);
            nmod_mpoly_mul(monomial_basis[c], monomial_basis[c], _var, ctx);
        }
        nmod_mpoly_mul(monomials_to_reduce[c], monomial_basis[c], xdegx, ctx);
    }
    nmod_mpoly_clear(_var, ctx);

    clock_t start, end_matrix, end_det;
    double cpu_time_used;

    start = clock();

    zz_p::init(MODULUS);

    Mat<zz_pX> M;
    M.SetDims(dimension_matrix, dimension_matrix);
    SetCoeff(M, degx, ident_mat_zz_p(dimension_matrix));

    clock_t column_start, column_end, normal_form_end;
    int var_index_read = 3;
    int var_index_write = 3;
    char filename_mult[40];
    snprintf(filename_mult, 40, "mult_mat_reduction_write%d.txt", nvars);
    FILE *flint_reduction_write = fopen(filename_mult, "w");
    #if !(COMPUTE_HARD_COLUMN)
    char filename_mult_read[40];
    snprintf(filename_mult_read, 40, "mult_mat_reduction_read.txt", nvars);
    nmod_mpoly_struct **HARD_COLUMNS = read_from_file_custom(filename_mult_read, ctx);
    int line = 0;
    #endif

    // fill in the matrix
    for (int j = 0; j < dimension_matrix; ++j)
    {
        column_start = clock();
        printf("computing column %d\n", j);
        
        /* We will compute the normal form of monomials_to_reduce[j], but cleverly... 
           we first compute the degree in each variable, assuming this is written on 2 words */

        mp_limb_signed_t degrees[nvars];
        nmod_mpoly_degrees_si(degrees, monomials_to_reduce[j], ctx);

        mp_limb_t index[nvars-1], sum_index = 0, var_index;
        for (int i = 0; i < nvars-1; ++i)
            index[i] = degrees[i] / weights[i];
        int multiply = 0;
        for (int i = 0; i < nvars-1; ++i)
        {
            if (degrees[i] != 0)
            {
                multiply = 1;
                var_index = i;
                index[i] -= 1;
                break;
            }
        }
        int fact = 1;
        for (int i = 0; i < nvars-1; ++i)
        {
            sum_index += fact*index[i];
            fact *= degs[i];
        }

        if (multiply)
        {
            nmod_mpoly_mul(monomials_to_reduce[j], monomials_to_reduce[sum_index], variables[var_index], ctx);
            puts("");
        }
        
        if(var_index >= var_index_read)
        {
            #if not COMPUTE_HARD_COLUMN
            monomials_to_reduce[j][0] = *HARD_COLUMNS[line];
            line+=1;
            #else
            
            normal_form(monomials_to_reduce[j], F, nvars, ctx);
            #endif
        }

        else
            normal_form(monomials_to_reduce[j], F, nvars, ctx);

        normal_form_end = clock();
        cpu_time_used = ((double) (normal_form_end - column_start)) / CLOCKS_PER_SEC;

        // nmod_mpoly_print_pretty(monomials_to_reduce[j], __x, ctx);
        printf("TIME FOR NORMAL FORM: %f s\n", cpu_time_used);
        fflush(stdout);
        mp_limb_t reduced_degs[nvars];
        mp_limb_t ind;
        mp_limb_t N = mpoly_words_per_exp(monomials_to_reduce[j]->bits, ctx->minfo);
        
        if(var_index>=var_index_write)
        {
            nmod_mpoly_fprint_pretty(flint_reduction_write, monomials_to_reduce[j], (const char**)__x, ctx);
            fprintf(flint_reduction_write,"\n");
            fflush(flint_reduction_write);
        }

        for (int i = 0; i < monomials_to_reduce[j]->length; ++i)
        {
            mpoly_unpack_vec_ui(reduced_degs, monomials_to_reduce[j]->exps + N*i, monomials_to_reduce[j]->bits, nvars, 2);
            ind = reduced_degs[0] / weights[nvars-1];
            int fact = 1, tmp = 0;
            for (int k = 0; k < nvars-1; ++k)
            {
                tmp += fact * reduced_degs[nvars-1-k] / weights[k];
                fact *= degs[k];
            }
            SetCoeff(M[tmp][j], ind, -zz_p(monomials_to_reduce[j]->coeffs[i]));
        }
    }

    fclose(flint_reduction_write);

    end_matrix = clock();
    cpu_time_used = ((double) (end_matrix - start)) / CLOCKS_PER_SEC;
    printf("MATRIX GENERATION TIME: %f\n", cpu_time_used);
    printf("Now computing determinant\n");
    fflush(stdout);
    // now the matrix has been computed
    zz_pX det;
    std::cout << "[det success] = " << determinant_generic_knowing_degree(det, M, degx* dimension_matrix) << std::endl;
    
    /* This function may be used instead of 'determinant_generic_knowing_degree' when genericity conditions are not met */
    // determinant_via_evaluation_geometric(det, M);

    zz_pX field_equation, _X;
    SetCoeff(field_equation, 1, zz_p(1));
    SetCoeff(_X, 1, zz_p(1));
    PowerMod(field_equation, field_equation, MODULUS, det);
    field_equation -= _X;


    det = GCD(det, field_equation);

    std::cout << det << std::endl;

    vec_zz_p roots = FindRoots(det);

    end_det = clock();
    cpu_time_used = ((double) (end_det - end_matrix)) / CLOCKS_PER_SEC;
    printf("DETERMINANT COMPUTATION TIME: %f\n", cpu_time_used);

    std::cout << "[SOLUTION]: " << roots << std::endl;

    /* Now write the solutions into a file */
    char filename[100];
    snprintf(filename, 100, "solutions_%d.txt", nvars);
    FILE *solutions = fopen(filename, "w");
    fprintf(solutions, "%ld\n", roots.length());
    for (int i = 0; i < roots.length(); ++i)
        fprintf(solutions, "%ld\n", conv<unsigned long, zz_p>(roots[i]));
    
    fclose(solutions);

    nmod_mpoly_clear(xdegx, ctx);

    for (int i = 0; i < dimension_matrix; ++i)
    {
        nmod_mpoly_clear(monomials_to_reduce[i], ctx);
        nmod_mpoly_clear(monomial_basis[i], ctx);
    }
        
    free(monomials_to_reduce);
    free(monomial_basis);

    // FREE_SYSTEM(F, ctx);
    for (int i = 0; i < nvars; ++i)
    {
        nmod_mpoly_clear(F[i], ctx);
        free(F[i]);
    }
    free(F);
    return 0;
}



void decomp_multi_index(slong *indices, slong *degs, int c, int nvars)
{
    for (int i = 0; i < nvars; ++i)
    {
        printf("degs[%d] = %ld\n", i, degs[i]);
        indices[i] = c % degs[i];
        c /= degs[i]; // DE
    }
}