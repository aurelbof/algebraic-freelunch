from time import time
import os
from tqdm import tqdm
import random
set_random_seed(64)


class Anemoi_Model:
    def __init__(self, p, alpha, l, n_r, g, round_constants=None, alpha_second=2):
        self.p = p
        self.F = GF(p)
        self.alpha = alpha
        self.alphainv = inverse_mod(alpha,p-1)
        self.l = l
        self.n_r = n_r
        self.g = F(g)
        self.ginv = 1/F(g)
        self.pol_system = []


        self.round_constants = round_constants
        if self.round_constants == None:
            self.round_constants = [self.F.random_element() for _ in range(2*l*n_r)]
        self.I = None
        self.alpha_second = alpha_second   # Exponent in H

        if l == 1:
            weighted_order = TermOrder('wdegrevlex',tuple(2 * 5**i * 3**(n_r - 1 - i) for i in range(n_r-1, -1, -1)) + (3**n_r,))   # comes from the formula wt(z_i^2) = wt(Q_i^3), and must be integers
            print(weighted_order)
        #TODO: l > 1 has a slightly different ordering
        

        print(weighted_order)
        # weighted_order=TermOrder('wdeglex',(18,30,50,27))
        # weighted_order='lex'

        self.R = PolynomialRing(self.F, n_r+1,
                                names=''.join(f'z{i},' for i in range(n_r-1, -1, -1)) + 'z' + str(n_r),    # z_n_r replaces x
                                order=weighted_order)
        self.z = list(self.R.gens_dict()[f'z{i}'] for i in range(n_r))
        self.x = self.R.gens_dict()['z' + str(n_r)]  # Otherwise x doesn't have the right parent
        #assert((self.x)**2 < (self.z[0])**alpha)

        if l == 1:
            self.Mx = matrix(F, 1, 1, [[1]])
            self.My = matrix(F, 1, 1, [[1]])

        #TODO: M if l is not 1


    def A(self, state, r):
        ret = deepcopy(state)
        for i in range(2* self.l):
            ret[i] += self.round_constants[2*(self.l)*r + i]
        return ret

    def M(self, state):
        ret = deepcopy(state)
        ret[: self.l] = self.Mx * ret[: self.l]
        ret[self.l :] = self.My * ret[self.l :]
        for i in range(self.l):
            ret[i+ self.l] += ret[i]
            ret[i] += ret[i+ self.l]
        return ret

    def sanitize(self, P):
        """Returns P reduced by all current polynomials in pol_system."""
        # return self.R(str(magma.NormalForm(P, self.pol_system)))
        # if self.QR == None:
        #     return P
        # return self.QR.lift(self.QR(P))
        if self.I == None: return P
        else: return P.reduce(self.I)
        
    def H(self, P, Q, r, model=False):
        PP = P - self.g * Q**self.alpha_second - self.ginv
        
        if model:
            #sanitize PP, as it may have zi's of high degree
            PP = self.sanitize(PP)
            self.pol_system.append(self.z[r]**(self.alpha) - PP)
            self.I = ideal(self.pol_system)
            QQ = Q - self.z[r]
        else:
            QQ = Q - PP**self.alphainv

        PP += g* QQ**self.alpha_second

        if model:
            #sanitize PP, as it may have zi's of high degree
            PP = self.sanitize(PP)

        return PP,QQ


    def S(self, state, r, model=False):
        ret = deepcopy(state)
        for i in range(self.l):
            ret[i],ret[i+ self.l] = self.H(ret[i], ret[i+ self.l], r, model=model)
        return ret

    def compute_model(self, homogeinize=True):
        start = time()
        if self.l == 1:
            state = vector(self.R, [self.x, 0]) 
        
        #TODO: cases where l > 1

        for r in tqdm(range(self.n_r)):
            state = self.A(state, r)
            state = self.M(state)
            state = self.S(state, r, model=True)
        state = self.A(state, r)
        state = self.M(state)
        end = time()
        print("time to compute polynomials before reduction:", end - start)
        filetime = f"anemoi/anemoi_gen_time_{self.n_r}.txt"
        if os.path.exists(filetime):
            os.remove(filetime)
        with open(filetime, 'x') as ft:
            ft.write(str(end - start) + '\n')
        # #CICO problem (0, *, ..., *)

        # start = time()
        # last_pol = state[0]
        # #Multiply last_pol by the right monomial in z_i and reduce to get a cool system
        # k_i = 1
        # for i in range(self.n_r):
        #     last_pol *= self.z[n_r - 1 - i]**((-k_i)%3)
        #     k_i = k_i + 2*ceil(k_i/3)
        # last_pol = self.sanitize(last_pol)

        # print("bad reduction time:", time() - start)


        # start = time()
        # last_pol = state[0]
        # #Multiply last_pol by the right monomial in z_i and reduce to get a cool system
        # k_i = 1
        # for i in range(self.n_r):
        #     last_pol *= self.z[n_r - 1 - i]**((-k_i)%3)
        #     last_pol = self.sanitize(last_pol)

        #     k_i = k_i + 2*ceil(k_i/3)
        # print("okay reduction time:", time() - start)

        last_pol = state[0]
        if homogeinize:
            start = time()
            #Multiply last_pol by the right monomial in z_i and reduce to get a cool system
            k_i = 1
            for i in tqdm(range(self.n_r)):
                for j in range((-k_i)%3):
                    last_pol *= self.z[self.n_r - 1 - i]
                    last_pol = self.sanitize(last_pol)
                k_i = k_i + 2*ceil(k_i/3)
            end = time()
            print("good reduction time:", end - start)
            filetime = f"anemoi/anemoi_sage_gstar_time_{self.n_r}.txt"
            if os.path.exists(filetime):
                os.remove(filetime)
            with open(filetime, 'x') as ft:
                ft.write(str(end - start) + '\n')
            last_pol = last_pol / (last_pol.lc())

        # start = time()
        # last_pol = state[0]
        # #Multiply last_pol by the right monomial in z_i and reduce to get a cool system
        # k_i = 1
        # k = []
        # for i in range(self.n_r):
        #     k.append((-k_i)%3)
        #     k_i += 2*ceil(k_i/3)
        # for i in tqdm(range(self.n_r)):
        #     for j in range(k[n_r-1-i]):
        #         last_pol *= self.z[i]
        #         last_pol = self.sanitize(last_pol)
        #     k_i = k_i + 2*ceil(k_i/3)
        # print("start with z0 reduction time:", time() - start)

        self.pol_system.append(last_pol)
        return self.pol_system

    def compute(self, state, number_of_rounds=None):
        if number_of_rounds == None:
            number_of_rounds = self.n_r
        for r in range(number_of_rounds):
            state = self.A(state, r)
            state = self.M(state)
            state = self.S(state, r)
        state = self.A(state, r)
        state = self.M(state)

        return state


def pol_system_to_file(pol_system, filename='anemoi_system.txt'):
    varnames = pol_system[0].parent().variable_names()
    nvars = len(varnames)
    varnames_reordered = [varnames[i] for i in range(nvars-2,-1,-1)] + [varnames[-1]]
    vars = pol_system[0].parent().gens()
    vars_reordered = [vars[i] for i in range(nvars-2,-1,-1)] + [vars[-1]]
    weights = pol_system[0].parent().term_order().weights()
    weights_reordered = [weights[i] for i in range(nvars-2,-1,-1)] + [weights[-1]]

    assert(pol_system[0].parent().term_order().is_weighted_degree_order())
    if os.path.exists('anemoi/' + filename):
        os.remove('anemoi/' + filename)
    with open('anemoi/' + filename, 'x') as f:
        f.write(str(nvars) + '\n')
        for i in range(nvars):
            f.write(varnames_reordered[i] + '\n')
        for i in range(nvars):
            f.write(str(weights_reordered[i]) + '\n')
        for p in pol_system:
            f.write(str(p).replace(' ', '') + '\n')

def pol_system_to_magma(pol_system, filename='anemoi.mag'):
    p = pol_system[0].parent().characteristic()
    varnames = pol_system[0].parent().variable_names()
    nvars = len(varnames)
    varnames_reordered = [varnames[i] for i in range(nvars-2,-1,-1)] + [varnames[-1]]
    vars = pol_system[0].parent().gens()
    vars_reordered = [vars[i] for i in range(nvars-2,-1,-1)] + [vars[-1]]
    weights = pol_system[0].parent().term_order().weights()
    weights_reordered = [weights[i] for i in range(nvars-2,-1,-1)] + [weights[-1]]

    if os.path.exists('anemoi/' + filename):
        os.remove('anemoi/' + filename)
    with open('anemoi/' + filename, 'x') as f:
        f.write('p := ' + str(p) + ';\n')
        f.write('R<' + str(varnames_reordered)[1:-1] + '>' + ':= PolynomialRing(FiniteField(p), ' + str(weights_reordered) + ');\n')
        f.write('eqns := ' + str(pol_system) + ';\n')
        f.write('n_r := ' + str(nvars-1) + ';\n')
        f.write('last_pol := eqns[n_r+1];\n')
        f.write('t := Realtime();\n')
        f.write('k_i := 1;\n')
        f.write('for i in [0..n_r-1] do\n')
        f.write('\tfor j in [0..(-k_i mod 3)-1] do\n')
        f.write('\t\tlast_pol := NormalForm(R.(n_r - i) * last_pol, eqns[1..n_r]);\n')
        f.write('\tend for;\n')
        f.write('\tk_i := k_i + 2*Ceiling(k_i/3);\n')
        f.write('end for;\n')
        f.write('time_taken := Realtime() - t;\n')
        f.write('SetAutoColumns(false);\n')
        f.write('SetColumns(0);\n')
        f.write('"Time to compute g*: " cat Sprint(time_taken);\n')
        f.write('Write("anemoi_gstar_time_" cat Sprint(n_r) cat ".txt", time_taken);\n')
        f.write('Write("anemoi_bigp_" cat Sprint(n_r) cat ".txt", Normalize(last_pol));\n')
        f.write('quit;\n')

    



if __name__ == "__main__":
    alpha = 3

    p = Primes().next(2**20)
    while (p-1)%alpha == 0:
        p = Primes().next(p+1)

    p = 28407454060060787

    print(f"p = {p}")

    alpha_inv = inverse_mod(alpha, p) 

    F = GF(p)

    l = 1
    g = 7
    solve = False


    for n_r in range(1,7):
        round_constants = [F.random_element() for _ in range(2*l*n_r)]
        if os.path.exists(f'anemoi/round_constants_{n_r}'):
            os.remove(f'anemoi/round_constants_{n_r}')
        with open(f'anemoi/round_constants_{n_r}', 'x') as f:
            f.write(str(round_constants))

        # anemoi = Anemoi_Model(p, alpha, l, n_r, g, round_constants=round_constants)
        # pol_system = anemoi.compute_model(homogeinize=True)
        anemoi_noh = Anemoi_Model(p, alpha, l, n_r, g, round_constants=round_constants)
        pol_system_noh = anemoi_noh.compute_model(homogeinize=False)
    
    
    # print("The polynomials are:")
    # for i in range(len(pol_system)):
    #     print(i, '\n',pol_system[i], '\n')
    # print()

        print('saving to file...')
        # pol_system_to_file(pol_system, filename = 'anemoi_system_' + str(n_r) + '.txt')
        pol_system_to_magma(pol_system_noh, filename = 'anemoi_' + str(n_r) + '.mag')
        pol_system_to_file(pol_system_noh, filename = 'anemoi_noh_' + str(n_r) + '.txt')

    if solve:
        I = ideal(pol_system)
        Ih = ideal(pol_system_noh)
        print("I is a grobner basis:", I.basis_is_groebner())

        V = I.variety()
        # V = [s for s in V if all(s.values())]   # remove parasitic soutions, i.e. ones where some z_i = 0
        print("solutions:", V)
        if not V:
            print("oops")
        for s in V:
            x_sol = s[pol_system[0].parent().gens()[-1]]
            anemoi_input = vector(F, [x_sol, 0])
            print("Anemoi on", anemoi_input, ":", anemoi.compute(anemoi_input))
        
        V = Ih.variety()
        # V = [s for s in V if all(s.values())]   # remove parasitic soutions, i.e. ones where some z_i = 0
        print("solutions:", V)
        if not V:
            print("oops")
        for s in V:
            x_sol = s[pol_system_noh[0].parent().gens()[-1]]
            anemoi_input = vector(F, [x_sol, 0])
            print("Anemoi on", anemoi_input, ":", anemoi.compute(anemoi_input))


        # print("\ntesting alt system:")
        # print(pol_system[:-1])
        # Ialt = ideal(pol_system[:-1] + [pol_system[0].parent().gens()[2]])
        # V = Ialt.variety()
        # V = [s for s in V if all(s.values())]   # remove parasitic soutions, i.e. ones where some z_i = 0
        # print("solutions:", V)
        # if not V:
        #     print("oops")
        # for s in V:
        #     x_sol = s[pol_system_noh[0].parent().gens()[-1]]
        #     anemoi_input = vector(F, [x_sol, 0])
        #     print("Anemoi on", anemoi_input, ":", anemoi.compute(anemoi_input))

    
    # print("F4 yields the following basis:")
    # G = ideal(pol_system).groebner_basis()
    # for i in range(len(G)):
    #     print(i, '\n',G[i], '\n')