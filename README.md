# How to find solutions to the CICO problem for GRIFFIN

These programs aim at find solutions for the CICO problem related to random instances of the GRIFFIN permutation with t=12 branches. The resolution proceeds as follows:

## STEP 0: installation requirements

Make sure that NTL and flint are both installed on the user's computer. The user must also download the 'Polynomial Matrix Library' available at 
`git@github.com:vneiger/pml.git`
and compile it wherever they want. Just make sure to know where the library is compiled. In the Makefile's located in generate\_system and solve\_system, change the value of the variable 'PATH\_TO\_PML' to the path of the folder ntl-extras on your system.

## STEP 1: generation of the polynomial equations

This step is implemented in the folder 'generate\_system'.
The user first executes
`sage input.sage <R>`
where 'R' is the number of rounds they wish to attack. The user may also change the value of the round constants by modifying line 15 which sets the random seed, but KEEP IN MIND that this implies a change in main.cpp in the last folder 'verify\_soltions' (this will be explained in STEP 3). The value of b at line 20 must not be change.

The user then lanches the command 'sage input.sage'. It will generate three files :
    - first\_two\_constants.txt which is located in the parent folder. Since 2 rounds of GRIFFIN are being bypassed by the attack.
    - constants%d.txt which is located in the current folder and stores the other round constants.
    - third\_state.txt which represents the internal state of the GRIFFIN permutation after 3 non-linear layers.

The program should compile using the command
`make`
Then the user launches the program with
`./main <nvars>`
Where nvars must be equal to 'R-2' (this is due to 3 rounds of Griffin being bypassed).
The program computes polynomial equations describing the CICO problem and write them in the file solve\_system/equations\_files/flint\_system{R-2}.txt.

## STEP 2: resolution of the equations

This step is implemented in the folder 'solve\_system'. The user specifies at line 30 of main.cpp the name of the file containing the equations. It should be equations\_files/flint\_system{R-2}.txt.
Compilation is done using
`make`
Execution is done using
`./main`
The program will generate a file named 'solutions\_{R-2}.txt' containing the solutions of the CICO problem (more exactly, the solutions for the last variable of the system of equations).

## STEP 3: verification of the solutions (optional)

The folder verify\_solutions provide a C++ implementation of GRIFFIN. To check the solutions, first run the file 'merge\_files.py'. This file takes nvars as an argument (reminder: nvars = R-2). this file generates a last file named 'solutions\_{R}\_rounds.txt'. In griffin.hpp, the user must hard-code the value of R at line 12. In main.cpp, the user must specify the name of the last file that was generated, i.e 'solutions\_{R}\_rounds.txt'.
Compilation is done using 'make', and execution using './main'. 

If the user changed the value of the round constants using 'input.sage' in 'generate\_system', the the values of 'line' (line 12) and 'consts' (line 13) must be modified accordingly: line corresponds to the coefficients in 'X' of the list of polynomials returned by 'input.sage', and 'consts' corresponds to the constant terms. For example, with the random seed given by the authors, 'input.sage' should return (among other things):

"""
Interesting line
0
6083092791329191\*X + 10262927513769948
73828876967109\*X + 26777978289746985
20242522495179319\*X + 25274272867876475
16618433409638936\*X + 26983928335905137
9524447280276718\*X + 26599315612406641
18608954680115099\*X + 27132545961972309
28215352360033571\*X + 17817841680635954
256760581109528\*X + 8525389239170846
4935843223313640\*X + 23893464510205431
2994222604449285\*X + 15186654911238877
23178583851352337\*X + 6203179491270662
"""

Therefore, lines 12 & 13 of main.cpp in verify\_solutions must be
    unsigned long long line[t] = {0, 6083092791329191, 73828876967109, 20242522495179319, 16618433409638936, 9524447280276718, 18608954680115099, 28215352360033571, 256760581109528, 4935843223313640, 2994222604449285, 23178583851352337};
    unsigned long long consts[t] = {0, 10262927513769948, 26777978289746985, 25274272867876475, 26983928335905137, 26599315612406641, 27132545961972309, 17817841680635954, 8525389239170846, 23893464510205431, 15186654911238877, 6203179491270662};
    
## Arion and Anemoi

Other files are provided in the folders 'Arion' and 'Anemoi' and implement the same kind of cryptanalysis against the two primitives.
