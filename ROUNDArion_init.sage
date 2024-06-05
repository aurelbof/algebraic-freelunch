import time
t1 = time.time() 
load("ArionHash.sage")
load("ROUNDArionHash.sage")
field = GF(28407454060060787)
branches = 3
rounds = 2
capacity = 2
d_1 = 3
d_2 = 121
matrix_type = 1
constants_g = [[2, 10], [9, 10], [4, 2], [3, 9]]
constants_h = [2, 6, 3, 9]
constants_aff = [[2, 6, 3], [7, 1, 5]]
arion_hash = ArionHash(field=field,
                             branches=branches,
                             rounds=rounds,
                             capacity=capacity,
                             d_1=d_1,
                             d_2=d_2,
                             constants_g=constants_g,
                             constants_h=constants_h,
                             constants_aff=constants_aff)
hash_val = 0
W = [d_2]+[3*d_1*(2**(branches-1)*(d_1+1)-d_1)**(i-1) + 1 for i in range(1,rounds)]
print("weight: " +str(W))

pols = generate_ArionHash_polynomials(field=field,
                                             branches=branches,
                                             rounds=rounds,
                                             capacity=capacity,
                                             d_1=d_1,
                                             d_2=d_2,
                                             constants_g=constants_g,
                                             constants_h=constants_h,
                                             constants_aff=constants_aff,
                                             hash_val=hash_val,
                                             order=TermOrder('wdeglex',W))
def generate_variables(branches, rate, rounds):
    variables = []
    for i in range(0, rate):
        variables.append("x")
    variables.append("z_" + str(1))
    for i in range(0, rounds - 2):
        variables.append("z_" + str(i + 2))
    return variables
variables=generate_variables(branches, branches-capacity, rounds)
P = PolynomialRing(field, variables)
file = open('Expeb.txt','w')
file.write(str(pols))
file.close()
t2 = time.time() 
print(t2-t1)
#print(pols)
#normal=P.ideal(pols).normal_basis()
#print(normal)
