#include <iostream>
#include <fstream>
#include <NTL/ZZ_p.h>
#include <NTL/vec_vec_ZZ_p.h>

#include "griffin.hpp"

using namespace std;

int main()
{
    unsigned long long line[t] = {0, 6083092791329191, 73828876967109, 20242522495179319, 16618433409638936, 9524447280276718, 18608954680115099, 28215352360033571, 256760581109528, 4935843223313640, 2994222604449285, 23178583851352337};
    unsigned long long consts[t] = {0, 10262927513769948, 26777978289746985, 25274272867876475, 26983928335905137, 26599315612406641, 27132545961972309, 17817841680635954, 8525389239170846, 23893464510205431, 15186654911238877, 6203179491270662};
    ZZ p(PRIME);
    ZZ_p::init(p);

    vec_vec_ZZ_p round_constants(INIT_SIZE, R);
    for (int i = 0; i < R; ++i)
        round_constants[i] = vec_ZZ_p(INIT_SIZE, t);

    std::ifstream solutions("solutions_5_rounds.txt");
    for (int i = 0; i < R; ++i)
        for (int j = 0; j < t; ++j)
            solutions >> round_constants[i][j];
    int nroots;
    solutions >> nroots;

    std::cout << "nroots = " << nroots << std::endl;

    vec_ZZ_p alphas(INIT_SIZE, t), betas(INIT_SIZE, t), gammas(INIT_SIZE, t);
    for (int i = 0; i < t; ++i)
    {
        alphas[i] = (i+1)*4;
        betas[i] = (i+1)*(i+1)*7;
        gammas[i] = i+1;
    }

    Griffin<ZZ_p> G(alphas, betas, gammas, ZZ(3));
    G.set_round_constants(round_constants);

    std::cout << G << std::endl;

    vec_ZZ_p input(INIT_SIZE, t);
    std::cout << "Beginning verif" << std::endl;
    ZZ_p x;
    for (int i = 0; i < nroots; ++i)
    {
        solutions >> x;
        std::cout << "testing " << x << std::endl;
        for (int j = 0; j < t; ++j)
            input[j] = line[j] * x + consts[j];
        vec_ZZ_p res = G(input);
        std::cout << "G(" << input << ") = " << res << std::endl;
    }
    return 0;
}