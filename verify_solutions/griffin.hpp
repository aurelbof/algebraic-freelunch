#ifndef GRIFFIN_HPP
#define GRIFFIN_HPP

#include <fstream>

#include <NTL/ZZ_p.h>
#include <NTL/mat_ZZ_pE.h>

using namespace NTL;

// the number of rounds/branches must be hard-coded 
#define R 5
#define t 12

#define PRIME 28407454060060787

#define SETNULL(M) for(int i = 0; i < t; ++i) \
                        for (int j = 0; j < t; ++j) \
                            M[i][j] = 0;


/**
 * @brief Implementation of the griffin permutation
 * 
 * @tparam T may be ZZ_p or ZZ_pE (i.e ZZ_pX mod some well chosen polynomial)
 */

template<class T>
class Griffin
{
    public:

        /**
         * @brief default constructor
         * d is chosen as the least odd prime s.t gcd(d, p-1) = 1
         * alpha and beta are chose randomly such that alpha**2 - 4*beta is not a square mod p
         */
        Griffin()
        {
            // ZZ_p alpha(1), beta(0);
            this->alpha = vec_ZZ_p(INIT_SIZE, t - 2);
            this->beta  = vec_ZZ_p(INIT_SIZE, t - 2);

            ZZ_p alpha2(1), beta2(0);
            ZZ_p disc = alpha2*alpha2 - 4*beta2;
            
            // first choose alpha2, beta2 s.t disc is not a square mod p
            while (Jacobi(disc.LoopHole(), ZZ(PRIME)) >= 0)
            {
                random(alpha2);
                random(beta2);
                disc = alpha2*alpha2 - 4*beta2;
            }

            for (int i = 0; i < t - 2; ++i)
            {
                this->alpha[i] = (i+1) * alpha2;
                this->beta[i]  = (i+1)*(i+1)*beta2;
            }
            
            this->gamma = vec_ZZ_p(INIT_SIZE, t - 2);
            for (int i = 0; i < t - 2; ++i)
            {
                this->gamma[i] = ZZ_p(0);
                while (this->gamma[i] == 0)
                    random(this->gamma[i]);
            }

            ZZ d(3);

            while (GCD(d, ZZ(PRIME) - 1) != 1)
                d = NextPrime(d+1);
            
            this->d = d;

            this->inv_d = InvMod(d, ZZ(PRIME) - 1);

            // random generation of round constants
            this->round_constants = Vec<Vec<T> >(INIT_SIZE, R);
            for (int i = 0; i < R-1; ++i)
                this->round_constants[i] = conv<Vec<T>, vec_ZZ_p>(random_vec_ZZ_p(t));

            // the last one is set to 0
            this->round_constants[R-1] = Vec<T>(INIT_SIZE, t);
            for (int i = 0; i < t; ++i)
                this->round_constants[R-1][i] = 0;

            // finally, build the MDS matrix
            build_mds();
        }

        /**
         * @brief basic constructor
         * round constants are randomly generated
         */
        Griffin(const vec_ZZ_p& alpha, const vec_ZZ_p& beta, const vec_ZZ_p& gamma, const ZZ& d)
        {
            this->alpha = alpha;
            this->beta  = beta;
            this->gamma = gamma;
            this->d     = d;

            this->inv_d = InvMod(d, ZZ(PRIME) - 1);

            // generation of the round constants
            this->round_constants = Vec<Vec<T> >(INIT_SIZE, R);
            for (int i = 0; i < R-1; ++i)
                this->round_constants[i] = conv<Vec<T>, vec_ZZ_p>(random_vec_ZZ_p(t));

            this->round_constants[R-1] = Vec<T>(INIT_SIZE, t);
            for (int i = 0; i < t; ++i)
                this->round_constants[R-1][i] = 0;

            build_mds();
        }

        /**
         * @brief build the MDS matrix
         * if t = 3, then M = M3 = [[2,1,1],[1,2,1],[1,1,2]]
         * if t = 4 then M = M4 = [[5,7,1,3],[4,6,1,1],[1,3,5,7],[1,1,4,6]]
         * if t = 4*t' >= 8 then M = [[2*M4  M4  ...  M4], [M4  2*M4  M4  ...  M4], ..., [M4  M4  ...  M4  2*M4]]
         */
        void build_mds()
        {
            this->MDS = Mat<T>(INIT_SIZE, t, t);

            if (t == 3)
            {
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j)
                        this->MDS[i][j] = i == j ? 2 : 1;
            }

            else
            {   
                // first build the matrix M'
                Mat<T> Mp(INIT_SIZE, t, t);
                SETNULL(Mp);
                int tp = t/4;

                int raw_M4[16] = {5,7,1,3,4,6,1,1,1,3,5,7,1,1,4,6};
                for (int b = 0; b < tp; ++b)
                    for (int i = 0; i < 4; ++i)
                        for (int j = 0; j < 4; ++j)
                            Mp[4*b + i][4*b + j] = raw_M4[4*i+j];

                if (t == 4)
                    this->MDS = Mp;

                else // then t = 4*t' >= 8
                {
                    Mat<T> Mpp(INIT_SIZE, t, t);

                    // build circ([2*I, I, ..., I])
                    for (int b1 = 0; b1 < tp; ++b1)
                        for (int b2 = 0; b2 < tp; ++b2)
                            for (int i = 0; i < 4; ++i)
                                for (int j = 0; j < 4; ++j)
                                    Mpp[4*b1+i][4*b2+j] = b1 == b2 ? (i == j ? 2 : 0) : (i == j ? 1 : 0);
                    
                    this->MDS = Mp * Mpp;
                }   
            }
        }

        void set_round_constants(const Vec<Vec<T> >& round_constants)
        {
            this->round_constants = round_constants;
        }

        void set_round_constant(const Vec<T>& round_constant, int i)
        {
            this->round_constants[i] = round_constant;
        }

        /**************************************
        *             PERMUTATION             *
        **************************************/

        // state = MDS * state + c[i]
        void linear_layer(Vec<T>& state, int i) const
        {
            mul(state, MDS, state);
            add(state, state, this->round_constants[i]);
        }


        void s_box(Vec<T>& state) const
        {
            Vec<T> x(state);
            state[0] = power(state[0], this->inv_d);
            state[1] = power(state[1], this->d);
            state[2] *= G_i(state[0], state[1], T(0), 0);
            for (int i = 3; i < t; ++i)
                state[i] *= G_i(state[0], state[1], x[i-1], i-2);
        }

        Vec<T> operator()(const Vec<T>& state) const
        {
            Vec<T> res(state);
            mul(res, this->MDS, state);
            for (int i = 0; i < R; ++i)
            {
                s_box(res);
                linear_layer(res, i);
            }
            return res;
        }

        Vec<T> get_round_constant(int i) const
        {
            return this->round_constants[i];
        }

        Vec<Vec<T> > get_round_constants() const
        {
            return this->round_constants;
        }

        ZZ get_d() const
        {
            return this->d;
        }

        ZZ get_inv_d() const
        {
            return this->inv_d;
        }

        vec_ZZ_p get_alpha() const
        {
            return this->alpha;
        }

        vec_ZZ_p get_beta() const
        {
            return this->beta;
        }

        vec_ZZ_p get_gamma() const
        {
            return this->gamma;
        }

        Mat<T> get_MDS() const
        {
            return this->MDS;
        }

        ZZ_pX get_modulus()
        {
            return this->__modulus;
        }

    private:

        T G_i(const T& x, const T& y, const T& z, int i) const
        {
            return power(gamma[i]*x + y + z, 2) + alpha[i] * (gamma[i]*x + y + z) + beta[i];
        }

        vec_ZZ_p alpha, beta, gamma;
        ZZ d, inv_d;
        Vec<Vec<T> > round_constants;
        Mat<T> MDS;
};


template<class T>
std::ostream& operator<<(std::ostream& stream, Griffin<T> g)
{
    stream << "exponant d = " << g.get_d() << std::endl;
    stream << "inverse = " << g.get_inv_d() << std::endl;
    stream << "alpha = " << g.get_alpha() << std::endl;
    stream << "beta = " << g.get_beta() << std::endl;
    stream << "gamma = " << g.get_gamma() << std::endl;
    stream << "round constants" << std::endl;
    for (int i = 0; i < R; ++i)
        stream << g.get_round_constant(i) << std::endl;

    stream << "Using matrix\n" << g.get_MDS() << std::endl;
    
    return stream;
}

#endif //GRIFFIN_HPP
