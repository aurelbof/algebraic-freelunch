#! /usr/bin/python3

import sys

def main(argv):
    if len(argv) != 2:
        print('Usage: ./merge_files.py <nvars>')
        return 1
    nvars = int(argv[1])
    
    first_two_constants = open('../first_two_constants.txt','r').read()
    constants_nvars = open(f'../generate_system/constants{nvars}.txt', 'r').read()
    solutions = open(f'../solve_system/solutions_{nvars}.txt').read()
    open(f'solutions_{nvars+2}_rounds.txt','w').write(first_two_constants + constants_nvars + solutions)

if __name__ == '__main__':
    main(sys.argv)