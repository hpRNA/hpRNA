################################################################################
##                                                                            ##
##  hamiltonian_path_generate.py                                      SCRIPT  ##
##  ------------------------------------------------------------------------  ##
##                                                                            ##
##  v. 08/10/2014                                                             ##
##                                                                            ##
##  James Geraets, University of York                                         ##
##  jg923@york.ac.uk                                                          ##
##                                                                            ##
##  Tool for generating a library of Hamiltonian paths corresponding to a     ##
##  particular geometry.                                                      ##
##                                                                            ##
################################################################################

### MODULE IMPORTS

import os
import sys
from string import maketrans

### GEOMETRY SPECIFICATION

# For ease, each of the 60 proteins can be labelled as a single character by
# using the following letters and numbers:
#   abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSUVTWXYZ01234567
# We have populated a connectivity map of neighbours. Two types of move are
# allowed. Either a move around a five_fold axis (either way) or across a
# two-fold axis.

# Proteins around a five_fold axes (vertex)
# C5 is "Clockwise 5": a move clockwise around a vertex
move_C5 = [('a', 't', 'N', 'u', 'd'),
           ('L', '6', 'R', 'O', 'M'),
           ('p', 'I', 'K', 's', 'r'),
           ('b', 'f', 'j', 'm', 'q'),
           ('c', 'w', 'y', 'g', 'e'),
           ('v', 'P', 'S', 'W', 'x'),
           ('Q', '2', 'Y', 'V', 'T'),
           ('G', '5', '3', '7', 'J'),
           ('l', 'F', 'H', 'o', 'n'),
           ('h', 'A', 'C', 'k', 'i'),
           ('z', 'X', 'U', '0', 'B'),
           ('D', '1', 'Z', '4', 'E')]

# Protein dimers
# DS is "Dimer Switch"
move_DS = [('a', 'b'),
           ('c', 'd'),
           ('e', 'f'),
           ('g', 'h'),
           ('i', 'j'),
           ('k', 'l'),
           ('m', 'n'),
           ('o', 'p'),
           ('q', 'r'),
           ('s', 't'),
           ('u', 'v'),
           ('w', 'x'),
           ('y', 'z'),
           ('A', 'B'),
           ('C', 'D'),
           ('E', 'F'),
           ('G', 'H'),
           ('I', 'J'),
           ('K', 'L'),
           ('M', 'N'),
           ('O', 'P'),
           ('Q', 'R'),
           ('S', 'T'),
           ('U', 'V'),
           ('W', 'X'),
           ('Y', 'Z'),
           ('0', '1'),
           ('2', '3'),
           ('4', '5'),
           ('6', '7')]

### FUNCTION DEFINITIONS

def clockwise(pos, move):
    for cyc in move:
        if pos in cyc:
            return cyc[(cyc.index(pos) + 1) % len(cyc)]

def anticlockwise(pos, move):
    for cyc in move:
        if pos in cyc:
            return cyc[(cyc.index(pos) - 1) % len(cyc)]

def find_path(iteration=None):
    '''
    Main function to find Hamiltonian paths. Lengthens each input path by every
    single move possible, lengthening by one. Can be restarted if the run is
    incomplete, if given the last fully successful iteration as an argument.
    
    For MS2, we only need paths of length 57, starting with a defined start of
    three proteins: 1,2,6 [a,b,f]. This reduces the number of calculated paths
    by symmetrical considerations. The full number can be recalculated
    afterwards. This is done in the hamiltonian_path_multiply.py script.
    
    '''
    
    if iteration == None:
        # For MS2, start at third iteration (symmetry considerations reduce
        # possible moves.
        iteration = 3
        # First three proteins are fixed (symmetry considerations).
        temp_list = ['abf']
        tempfile = open('paths/paths_%02i.txt' % (iteration,), 'w')
        tempfile.write('\n'.join(temp_list))
        tempfile.close()
    
    # Only paths of length 57 are required - special rules for MS2 ignoring
    # five-fold (short) moves around start/end vertex.
    # Not all paths output will be cyclic - further step required.
    # Some output paths will violate special start/end rules for MS2. These
    # also need removing (further step).
    while iteration < 58:
        infile = open('paths/paths_%02i.txt' % (iteration,), 'r')
        
        iteration += 1
        print iteration
        
        outfile = open('paths/paths_%02i.txt' % (iteration,), 'w')
        
        for i in infile:
            execute(i.strip(), outfile)
        
        infile.close()
        outfile.close()


def execute(i, outfile):
    '''
    Three moves. Extend each path by one unit, attempting each possible move.
    
    '''
    # Move 1
    np = m1(i[-1])
    if np not in i:
        outfile.write(i+np+'\n')
    # Move 2
    np = m2(i[-1])
    if np not in i:
        outfile.write(i+np+'\n')
    # Move 3
    np = m3(i[-1])
    if np not in i:
        outfile.write(i+np+'\n')

def m1(pos):
    '''
    Dimer switch is move 1.
    
    '''
    return clockwise(pos, move_DS)
    
def m2(pos):
    '''
    Clockwise around five-fold is move 2.
    
    '''
    return clockwise(pos, move_C5)

def m3(pos):
    '''
    Anticlockwise around five-fold is move 3.
    
    '''
    return anticlockwise(pos, move_C5)

### MAIN

if __name__ == '__main__':
    if not os.path.exists('paths'):
        os.makedirs('paths')    # Make a new directory  
    find_path()
    
## ENDS