################################################################################
##                                                                            ##
##  hamiltonian_path_instantiate.py                                   SCRIPT  ##
##  ------------------------------------------------------------------------  ##
##                                                                            ##
##  v. 08/10/2014                                                             ##
##                                                                            ##
##  James Geraets, University of York                                         ##
##  jg923@york.ac.uk                                                          ##
##                                                                            ##
##  Tool for generating instances of Hamiltonian paths satisfying viral       ##
##  geometries, starting and ending at specified positions.                   ##
##                                                                            ##
################################################################################

### MODULE IMPORTS

from math import sqrt
import os.path
import os
import pickle
import sys

### GEOMETRY SPECIFICATION

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


# Proteins where the paths may start. Constrained to be close to the maturation
# protein of MS2.

starting_proteins = [('a', 't', 'N', 'u', 'd'),
                     ('L', '6', 'R', 'O', 'M'),
                     ('p', 'I', 'K', 's', 'r'),
                     ('b', 'f', 'j', 'm', 'q'),
                     ('v', 'P', 'S', 'W', 'x'),
                     ('Q', '2', 'Y', 'V', 'T'),
                     ('G', '5', '3', '7', 'J'),
                     ('l', 'F', 'H', 'o', 'n')]

### FUNCTION DEFINITIONS

def score_long_gen():
    return dict((a, 0) for a in move_DS)
    
def score_short_gen():
    score_short = dict()
    for v in move_C5:
        for i in range(5):
            score_short[(v[i-1],v[i])] = 0
    return score_short

def clockwise(pos, move):
    '''
    Cycle through the moves clockwise.
    
    '''
    for cyc in move:
        if pos in cyc:
            return cyc[(cyc.index(pos) + 1) % len(cyc)]

def anticlockwise(pos, move):
    '''
    Cycle through the moves anticlockwise.
    
    '''
    for cyc in move:
        if pos in cyc:
            return cyc[(cyc.index(pos) - 1) % len(cyc)]
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

def hami(string_path, start_point, longdict, shortdict):
    position = start_point
    for i in string_path:
        if i == '1':
            temp = m1(position)
            if (position, temp) in longdict:
                longdict[(position, temp)] = 1
            else:
                longdict[(temp, position)] = 1
            position = temp
        elif i == '2':
            temp = m2(position)
            if (position, temp) in shortdict:
                shortdict[(position, temp)] = 1
            else:
                shortdict[(temp, position)] = 1
            position = temp
        elif i == '3':
            temp = m3(position)
            if (position, temp) in shortdict:
                shortdict[(position, temp)] = 1
            else:
                shortdict[(temp, position)] = 1
            position = temp

    return longdict, shortdict
    
def load_hamiltonian_paths(filename):
    '''
    Load in paths in 1,2,3 notation.
    
    '''
    filein = open(filename, 'r')
    # Set construct allows any duplicates to be removed
    paths = set()
    for line in filein:
        paths.add(line.strip())
    filein.close()
    return paths
    
def main(hampaths):
    
    realised_paths = list()
    consolidation_set = set()
    
    for vertex in starting_proteins:
        for start in vertex:
        
            for index, hampath in enumerate(hampaths):
                print index
                print hampath
                longdict, shortdict = hami(hampath, start, score_long_gen(), score_short_gen())
        
                realised_paths.append((hampath, start, longdict, shortdict))
                consolidation_set.add((''.join(str(stuff[1]) for stuff in sorted(longdict.items())), ''.join(str(stuff[1]) for stuff in sorted(shortdict.items()))))
    
    return realised_paths, consolidation_set

### MAIN

if __name__ == '__main__':
    assert len(sys.argv) == 2
    
    hpath_input_name, hpath_input_extension = os.path.splitext(sys.argv[1])
    hampaths = load_hamiltonian_paths(hpath_input_name + '.hpath')
    realised_paths, consolidation_set = main(hampaths)
    
    realiseddump = open(hpath_input_name + '.pick', 'w')
    pickle.dump(realised_paths, realiseddump)
    realiseddump.close()
    
    print '\n\n'
    print 'CHECK UNIQUE'
    print len(consolidation_set), len(realised_paths)
    #print consolidation_set

### ENDS