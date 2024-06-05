"""
    Polynomial model of arion_hash.
"""

load("ArionHash.sage")

def generate_variables(branches, rate, rounds):
    variables = []
    for i in range(0, rate):
        variables.append("x")
    variables.append("z_" + str(1))
    for i in range(0, rounds - 2):
        variables.append("z_" + str(i + 2))
    return variables

def evaluate_g_and_h(x_in, constants_g, constant_h):
    # Constant term
    out_g = constants_g[1]
    out_h = 0
    # Linear term
    out_g += x_in * constants_g[0]
    out_h += x_in * constant_h
    # Quadratic term
    x_temp = x_in**2
    out_g += x_temp
    out_h += x_temp
    return out_g, out_h

def GTDS(v_in, v_n_std, d_1, constants_g, constants_h):
    branches = len(v_in)
    v_out = copy(v_in)
    sigma = v_n_std + v_out[branches - 1]
    for i in range(branches - 2, -1, -1):
        v_out[i] = v_out[i]**d_1
        g_i, h_i = evaluate_g_and_h(sigma, constants_g[i], constants_h[i])
        v_out[i] *= g_i
        v_out[i] += h_i
        sigma += v_in[i] + v_out[i]
    return v_out

def affine_layer(v_in, matrix, round_constants):
    return matrix * v_in + round_constants

def round_function(v_in, v_n_std, d_1, matrix, constants_g, constants_h, constants_aff):
    v_out = GTDS(v_in, v_n_std, d_1, constants_g, constants_h)
    print("GTDS done")
    v_out = affine_layer(vector(v_out), matrix, constants_aff)
    print("affine layer done")
    return v_out
    
def W(branches, rounds):
    W = [d_2] + [3*d_1*(2**(branches-1)*(d_1+1)-d_1)**(i-1) + 1 for i in range(1,rounds)]

print(W)

def generate_ArionHash_polynomials(field=GF(0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001), # BLS12
                                   branches=3,
                                   rounds=6,
                                   capacity=2,
                                   d_1=None,
                                   d_2=None,
                                   constants_g=None,
                                   constants_h=None,
                                   constants_aff=None,
                                   hash_val=None,
                                   order=TermOrder('wdeglex',W),
                                   naive_model=False):
    
    arion_hash = ArionHash(field=field,
                           branches=branches,
                           rounds=rounds,
                           capacity=capacity,
                           d_1=d_1,
                           d_2=d_2,
                           constants_g=constants_g,
                           constants_h=constants_h,
                           constants_aff=constants_aff)
    print_plain = False
    if hash_val is None:
        plain = field.random_element()
        hash_val = arion_hash.hash(plain)
        print_plain = True
    
    if print_plain:
        print("Plain text: " + str(plain))
    print("Hash value: " + str(hash_val))
    print("Term order: " + str(order))
    
    variables = generate_variables(branches, arion_hash.rate, rounds)
    P = PolynomialRing(field, variables, order=order)
    variables = [P(var) for var in variables]
    polynomials = []
    pols=[]
    counter = arion_hash.rate
    
    matrix=arion_hash.matrix
    vari = [1]
    for i in range(1, branches-1):
        vari.append(var("c_" + str(i)))
    PP=PolynomialRing(field, variables+vari[1:])
    varii=[1]+[PP(var) for var in vari[1:]]+[0]
    new_vector = matrix * vector(varii)
    print(new_vector)
    new_vector=[new_vector[j] for j in range(2,branches)]
        
    #lineq=solve((new_vector[i]==0 for i in range(0,branches-2)), vari[1:])
    #lineq=solve((3*c_1+2==0), vari[1:])
    #print("lineq:"+str(lineq))
    
    #for i in range(1,branches-1):
    #    varii[i]=lineq[i-1].right()
    #varii=[1,P(-2/3),0]
    #varii=[1,P(-10/13),P(1/13),0]
    #varii=[1,P(-5/6),P(1/18),P(1/18),0]
    #varii=[1,P(-7/8),P(1/24),P(1/24),P(1/24),0]
    varii=[1,P(-12/13),P(15/585),P(15/585),P(1/39),P(1/39),P(1/39),0]
    current_state=[]
    for i in range(0,branches):
        current_state.append(varii[i]*P(x))  
    
    tmp = current_state[branches - 1]
    print("initial state" + str(current_state))
    if rounds > 1:
        current_state = affine_layer(vector(current_state), matrix, 0)
        print("after matrix"+str(current_state))
        polynomials.append(list(round_function(current_state,
                                               tmp,
                                               arion_hash.d_1,
                                               arion_hash.matrix,
                                               arion_hash.constants_g[0:branches - 1],
                                               arion_hash.constants_h[0:branches - 1],
                                               vector(field, arion_hash.constants_aff[0]))))
        #print("polynomials first round"+str(polynomials))
        for r in range(1, rounds):
            print("starting round "+str(r+1))
            current_state=vector(polynomials[0])
            tmp = current_state[branches - 1]
            current_state[branches - 1] = variables[counter]
            polynomials=[]
            polynomials.append(list(round_function(current_state,
                                                   tmp,
                                                   arion_hash.d_1,
                                                   arion_hash.matrix,
                                                   arion_hash.constants_g[r * (branches - 1):(r + 1) * (branches - 1)],
                                                   arion_hash.constants_h[r * (branches - 1):(r + 1) * (branches - 1)],
                                                   vector(field, arion_hash.constants_aff[r]))))
            pols.append(tmp - current_state[branches - 1]**arion_hash.d_2)
            counter+=1
        polynomials=polynomials[0]
        polynomials=polynomials[0]
        pols.append(polynomials)
    
    return pols
