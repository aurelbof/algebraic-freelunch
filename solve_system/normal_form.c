#include <stdlib.h>

#include "normal_form.h"

void normal_form(nmod_mpoly_t h, nmod_mpoly_struct **F, int N, nmod_mpoly_ctx_t ctx)
{
    // printf("debug: N = %ld\n", N);
    nmod_mpoly_struct **Q = (nmod_mpoly_struct **) malloc(N*sizeof(nmod_mpoly_struct*));
    for (int i = 0; i < N; ++i)
    {
        Q[i] = (nmod_mpoly_struct *) malloc(sizeof(nmod_mpoly_struct));
        nmod_mpoly_init(Q[i], ctx);
    }

    nmod_mpoly_t H;
    nmod_mpoly_init(H, ctx);
    nmod_mpoly_set(H, h, ctx);
    nmod_mpoly_divrem_ideal(Q, h, H, F, N, ctx);

    
    for (int i = 0; i < N; ++i)
    {
        nmod_mpoly_clear(Q[i], ctx);
        free(Q[i]);
    }
    free(Q);                         
    nmod_mpoly_clear(H, ctx);  
}