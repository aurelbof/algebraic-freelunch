# How to compile
To compile the program, just type the command line 'make' in a terminal.
Make sure to have NTL (and gmp) installed on your computer.

# How to change the number of rounds
Change the value of 'R' at line 12 in 'griffin.hpp' and the name of the file 'solutions' in main.cpp at line 21. Then recompile.

# Format of the files 'solutions_{5,6,7}.txt'
The program expects the file 'solutions{5,6,7}.txt' to have the following format:
constant[0]
constant[1]
...
constant[r-1]
number of solutions n
solution_1
...
solution_n

Each constant is actually a vector of dimension t=12.

# format of the output
The program outputs some data about the Griffin instance that is being used.
Below the line 'beginning verif', the program prints 'testing {value}' for each value provided in 'solutions_{5,6,7}.txt'. It then prints G([0 * * ... *]) = [* * * ... *] and you should check that the first coordinate of the output is 0.