################################################################################
##                                                                            ##
##  hamiltonian_path_constrain.py                                     SCRIPT  ##
##  ------------------------------------------------------------------------  ##
##                                                                            ##
##  v. 08/10/2014                                                             ##
##                                                                            ##
##  James Geraets, University of York                                         ##
##  jg923@york.ac.uk                                                          ##
##                                                                            ##
##  Tool for filtering instances of realised Hamiltonian paths on 3-d         ##
##  geometries, fulfilling constraints on routes.                             ##
##                                                                            ##
################################################################################

### MODULE IMPORTS

import sys
import os
import pickle
import string
import copy
from hamiltonian_path_draw import hami_draw

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

def permute_path(path):
    '''
    Return all geometries of a path.
    
    '''
    return [path,
            backwards_string(path),
            mirror_string(path),
            reverse_string(path)]

def backwards_string(path):
    '''
    Returns a path backwards, creating a mirror path.
    
    '''
    return path[::-1]

def mirror_string(path):
    '''
    Translates a path such that all instances of 2 become 3, and all instances
    of 3 become 2. Creates a mirror path.
    
    '''
    trans = string.maketrans('23','32')
    return path.translate(trans)
    
def reverse_string(path):
    '''
    Translates a path such that all instances of 2 become 3, and all instances
    of 3 become 2, and then backs up the path. Creates a mirror path, that is
    the equivalent of describing the same path in reverse, i.e. from 3' to 5'
    instead of 5' to 3'.
    
    '''
    return mirror_string(backwards_string(path))
        
def dump_hampaths(paths):
    dfilename, dfileextension = os.path.splitext(sys.argv[1])
    rfilename, rfileextension = os.path.splitext(sys.argv[2])
    pidump = open(dfilename + '.' + rfilename + '.pick', 'w')
    pickle.dump(paths, pidump)
    pidump.close()

def comparison(input_paths, output_paths):
    '''
    Print a little comparison before and after constraints
    
    '''
    dump_hampaths(output_paths)
    
    best_paths = count_best_paths(output_paths)
    str1 = '**   original:  ' + str(len(input_paths))
    str2 = '**   processed: ' + str(len(output_paths))
    str3 = '**   best:      ' + str(best_paths)
    str1 += (' '*(28-len(str1))) + '**'
    str2 += (' '*(28-len(str2))) + '**'
    str3 += (' '*(28-len(str3))) + '**'
    
    print '\n******************************'
    print '**                          **'
    print str1
    print str2
    print str3
    print '**                          **'
    print '******************************\n'

def main(pick_data, constraints_filename):
    
    # Load in the traced Hamiltonian paths in pickled format
    piload = open(pick_data, 'r')
    input_paths = pickle.load(piload)
    piload.close()
    lentot = len(input_paths)

    # Load in file with constrained edges
    # Sample constraint file provided
    constraints_file = open(constraints_filename)
    output_paths = list()
    constrain_long = dict()
    constrain_short = dict()
    for line in constraints_file:
        spline = line.strip().split()
        # Constraints have 4 elements
        assert len(spline) == 4
        
        # Long edge constraints
        if spline[0] == 'L':
            constrain_long[(spline[2], spline[3])] = int(spline[1])

        # Short edge constraints
        elif spline[0] == 'S':
            constrain_short[(spline[2], spline[3])] = int(spline[1])
    constraints_file.close()

    # Sift through paths, removing those that do not meet constraints
    for hampath, startname, longdict, shortdict in input_paths:
        if all(longdict[(n1, n2)] == e for (n1, n2), e in constrain_long.iteritems()):
            if all(shortdict[(n1, n2)] == e for (n1, n2), e in constrain_short.iteritems()):
                output_paths.append((hampath, startname, longdict, shortdict))
    output_paths = upshift(output_paths)
    if len(output_paths) < 20 and len(output_paths) > 0:
        # If few result paths, display and draw.
        pngname = '-'.join(pick_data.split('.')[:2]) + '-' + os.path.splitext(constraints_filename)[0]
        display_solution_paths(input_paths, output_paths[:5], pngname, constrain_long, constrain_short)
        cl, cs = inferred(output_paths, constrain_long, constrain_short)
        
    elif len(output_paths) > 0:
        # Otherwise, just list them.
        comparison(input_paths, output_paths)
        cl, cs = inferred(output_paths, constrain_long, constrain_short)
        
    else:
        # If no result paths, the constraints you have imposed mean none are feasible.
        print '\nNO SOLUTIONS'

def inferred(output_paths, constrain_long, constrain_short):
    '''
    Checks to see if constraints has mandated certain routes.
    Returns edges, long and short, that are always the same in the remaining 
    solutions. Either they are always used, or always not used.
    So these are inferences of the geometry and constraints set.
    
    '''
    check_hampath, check_startname, check_longdict, check_shortdict = output_paths[0]


    for hampath, startname, longdict, shortdict in output_paths:

        check_longdict = dict((k,v) for k, v in check_longdict.iteritems() if v == longdict[k])
        check_shortdict = dict((k,v) for k, v in check_shortdict.iteritems() if v == shortdict[k])

    check_longdict = dict((m, n) for m, n in check_longdict.items() if (m, n) not in constrain_long.items())
    check_shortdict = dict((m, n) for m, n in check_shortdict.items() if (m, n) not in constrain_short.items())

    if len(check_longdict) > 0 or len(check_shortdict) > 0:
        print '\n***\nINFERRED EDGES'
        print '\nLONG:  '
        for (n1, n2), e in sorted(check_longdict.iteritems()):
            print n1 + '\t' + n2 + '\t' + str(e)
        print '\nSHORT: '
        for (n1, n2), e in sorted(check_shortdict.iteritems()):
            print n1 + '\t' + n2 + '\t' + str(e)
        print '***\n'
    else:
        print '\nNO INFERRED EDGES\n'
    return check_longdict, check_shortdict
    
def display_solution_paths(input_paths, output_paths, pngname, constrain_long, constrain_short):
    
    comparison(input_paths, output_paths)
    
    print 'SOLUTION PATHS\n'
    for hampath, startname, longdict, shortdict in output_paths:
        draw_out_long = dict((k,v) for k,v in longdict.iteritems() if k not in constrain_long)
        draw_out_short = dict((k,v) for k,v in shortdict.iteritems() if k not in constrain_short)
        print startname, hampath
        hami_draw(constrain_long, constrain_short, draw_out_long, draw_out_short, pngname + '-' + hampath + '-' + startname)

def count_best_paths(output_paths):
    '''
    Count the remaining paths that correspond to the best Hamiltonian path.
    Either of the four orientations of it left by hamiltonian_path_generation.py
    
    '''
    best_out = []
    
    for best_p in permute_path('13133331333313312222133331222213333122122221312133331331'):
        best_out += filter(lambda x: x[0] == best_p, output_paths)

    return len(best_out)
    
def upshift(paths_in, instances=None):
    '''
    Moves paths corresponding to instances to top of the list, for easy finding
    manually when/if printed to file.
    
    '''
    if instances == None:
        instances = ['12213313133121333313133331331213122221221333312222133331',
                     '12122221221212222131221212213312222133121312212122133131',
                     '13313133331333313312131222212213333122221333312222133331',
                     '13133331333313312222133331222213333122122221312133331331']
    
    full_instances = []
    for instance in instances:
        full_instances.extend(permute_path(instance))
    
    paths_out_first = []
    paths_out_second = []
    for path in paths_in:
        if path[0] in full_instances:
            paths_out_first.append(path)
        else:
            paths_out_second.append(path)
    
    return paths_out_first + paths_out_second

### MAIN

if __name__ == '__main__':
    # Three arguments must be supplied
    # Usage: python refine.py [pickled data file] [constraints file]
    assert len(sys.argv) == 3
    pick_data = sys.argv[1]
    constraints_filename = sys.argv[2]
    main(pick_data, constraints_filename)
    
### ENDS