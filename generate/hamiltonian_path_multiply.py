################################################################################
##                                                                            ##
##  hamiltonian_path_multiply.py                                      SCRIPT  ##
##  ------------------------------------------------------------------------  ##
##                                                                            ##
##  v. 08/10/2014                                                             ##
##                                                                            ##
##  James Geraets, University of York                                         ##
##  jg923@york.ac.uk                                                          ##
##                                                                            ##
##  Tool for validating generated Hamiltonian paths, and generating           ##
##  symmetrically equivalent Hamiltonian paths to form a library. Cyclic      ##
##  paths are saved, that start and end from the same vertex. Paths are not   ##
##  permitted to visit start/end vertex during middle of path. The generated  ##
##  library can be used to interrogate viral RNA structure information.       ##
##  Paths are represented in 1,2,3 notation: moves between proteins are       ##
##  notated rather than the proteins themselves. The generalised notation     ##
##  is more effective for identifying symmetries.
##                                                                            ##
################################################################################

### MODULE IMPORTS

import string

### GEOMETRY SPECIFICATION

# The 60 proteins can be labelled as a single character:
#   abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSUVTWXYZ01234567

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

def wrapin(path, ind):
    '''
    Function to split a 'cyclic' path at a given index, recombining so the
    index point becomes the start of the path.
    
    '''
    return path[ind:] + path[:ind]

def wrapout(path, ind):
    '''
    Function to split a 'cyclic' path at a given index, recombining so the
    index point is removed from the path.
    
    '''
    return path[ind + 1:] + path[:ind]
    
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
    
def trim_path(path):
    '''
    Removes 2 and 3 numbers from the beginning and end of a path. This means
    the path that the path represents ignores the initial and final small
    moves. In our analysis, we are unable to determine the small moves, so
    the multiple options for starting/finishing the paths are redundant.
    
    '''
    return path[path.index('1'):path.rindex('1') + 1]
    
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
    
def save_hamiltonian_paths(paths, filename):
    '''
    Write out paths in 1,2,3 notation.
    
    '''
    fileout = open(filename, 'w')
    for path in paths:
        fileout.write(path + '\n')
    fileout.close()

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

def number_system(path):
    '''
    Converts a protein-notated path to 1,2,3 notation. 1 represents a dimer
    switch, 2 a clockwise move around a five-fold, and 3 an anticlockwise move.
    
    '''
    str = ""
    for i, pos in enumerate(path[:-1]):
        if m1(pos) == path[i + 1]:
            str += "1"
        elif m2(pos) == path[i + 1]:
            str += "2"
        elif m3(pos) == path[i + 1]:
            str += "3"
        else:
            # Raise an error if not neighbours.
            assert False
    return str
    
def load_protein_paths(filename):
    '''
    Loads in generated Hamiltonian paths with specified protein order. Converts
    to generalised 1,2,3 notation referring to types of moves rather than
    specific proteins. This allows paths to be compared and computed quicker.
    Also finds degenerate symmetrical paths by calculation. Sifts for circular
    paths only.
    
    '''
    paths = []
    infile = open(filename, 'r')
    for line in infile:
        paths.append(line.strip())
    infile.close()
    
    # We want only paths that start/end at same vertex as 'a'
    # included for completeness
    endings = set(['a', 't', 'N', 'u', 'd'])
    
    # Also, special rules - allowed us to only generate up to length 57 instead
    # of length 59: we do not allow Hamiltonian path to visit the start/end
    # vertex at any other point. We test to see if sets are disjoint.
    circular_paths = set()
    for path in paths:
        # Check end condition, not visit primary vertex condition
        if path[-1] in endings and endings.isdisjoint(path[1:-1]):
            numpath = number_system(path)
            # Add all symmetries at this stage for transparency
            circular_paths.add(numpath)
            circular_paths.add(backwards_string(numpath))
            circular_paths.add(mirror_string(numpath))
            circular_paths.add(reverse_string(numpath))
    return circular_paths

def main_primary():
    '''
    Use our generated protein-notated paths to form 
    '''
    protein_paths = load_protein_paths('paths_57.txt')
    
    print 'protein', len(protein_paths)
    
    # There are 264 unique paths.
    # However, as each path starts and ends at two different valid starting and
    # ending locations, we have a two-fold degeneracy, as each path can be
    # equally validly described in each direction. I.e., a path can be described
    # from its 3' end or 5' end in the algorithm equally. We thus have to remove 
    # half of the paths.
    
    final_paths = set()
    
    for path in sorted(list(protein_paths)):
        revpath = reverse_string(path)
        
        # If the reverse path is not in final set, add the path.
        if revpath not in final_paths:
            final_paths.add(path)
    
    # After removing this degeneracy, there are 132 paths.
    # Each of these paths can start at any packaging signal position nominated.
    # For our analysis, we constrained the start and end of the paths to be at
    # the eight vertices nearest the maturation protein of MS2.
    
    # This meant the total realisation of paths on the sphere was 132*5*8=5,280
    
    # Note that this includes the fact that in our analysis, 5' and 3'
    # directions are not distinguishable.
    
    print 'final', len(final_paths)
    save_hamiltonian_paths(final_paths, 'final_paths.hpath')

def main_secondary():
    '''
    Link to previous quoted Hamiltonian path libraries for MS2. Uses libraries
    (txt files) for cycles and pseudo cycles on the MS2 geometry to generate
    paths.
    
    '''
    # Load in the extant libraries of Hamiltonian paths:
    # Full cyclic Hamiltonian paths
    cycle_library = load_hamiltonian_paths('ham_cycles.txt')
    # So-called "pseudo-cycle" Hamiltonian paths
    pcycle_library = load_hamiltonian_paths('ham_pseudocycles.txt')
    
    # Permute the library paths, as they only contain one symmetry rotation of
    # each path.
    cycle_paths = set()
    pcycle_paths = set()
    
    # In a cycle, every time there is a full loop around a vertex, this could
    # be replaced by a start/end special vertex (for MS2).
    # This calculates the realisations of the cycle through the MS2 case,
    # leaving intact the short edges around the five-fold (barring a single
    # break each time). This equivalizes to the pseudo-cycles (of length 59).
    for c in cycle_library:
        d = c[:]
        for i in range(c.count('2222')):
            ind = d.index('2222')
            cycle_paths.add(wrapout(d, ind))
            cycle_paths.add(wrapout(d, ind + 1))
            cycle_paths.add(wrapout(d, ind + 2))
            cycle_paths.add(wrapout(d, ind + 3))
            
            d = wrapin(d, ind + 5)
            
        d = c[:]
        for i in range(c.count('3333')):
            ind = d.index('3333')
    
            cycle_paths.add(wrapout(d, ind))
            cycle_paths.add(wrapout(d, ind + 1))
            cycle_paths.add(wrapout(d, ind + 2))
            cycle_paths.add(wrapout(d, ind + 3))
            
            d = wrapin(d, ind + 5)
            
    print 'cycles', len(cycle_paths)
    
    # Pseudo cycles on the other hand, do not complete a full cycle.
    # They can only start and end at one vertex. However there are two
    # permutations for the short edges for each pseudo cycle.
    # In the text file, only 22... ...3 paths are listed.
    #                        3... ...22 paths are alternate.
    for p in pcycle_library:
        pcycle_paths.add(p)
        pcycle_paths.add('3' + p[2:-1] + '22')
    
    print 'pcycles', len(pcycle_paths)
    
    # Cycles and pseudo cycles combined.
    
    all_paths = set()
    
    for path in (cycle_paths | pcycle_paths):
        all_paths.add(path)
        all_paths.add(mirror_string(path))
        all_paths.add(backwards_string(path))
        all_paths.add(reverse_string(path))
        
    print 'total', len(all_paths)
    save_hamiltonian_paths(all_paths, 'all_paths.hpath')
    
    # For the supplied libraries, there are 824 total paths found.
    # Each of the 824 paths can start or end at any of 5 packaging signal
    # postions at any of 12 vertices. This leads to 824*5*12 = 49,440 possible
    # path realisations for MS2, quoted previously.
    
    # In our analysis, we cannot determine short edges. Thus short edges not
    # essential to determine the path, i.e. the beginning and ending moves, are
    # not needed for the analysis. These are trimmed.
    
    trimmed_paths = set()
    for path in list(all_paths):
        trimmed_paths.add(trim_path(path))
    
    print 'trimmed', len(trimmed_paths)
    save_hamiltonian_paths(trimmed_paths, 'trimmed_paths.hpath')
    
    # After trimming there are 264 paths.
    # However, as each path starts and ends at two different valid starting and
    # ending locations, we have a two-fold degeneracy, as each path can be
    # equally validly described in each direction. I.e., a path can be described
    # from its 3' end or 5' end in the algorithm equally. We thus have to remove 
    # half of the paths.
    
    individual_paths = set()
    
    for path in sorted(list(trimmed_paths)):
        revpath = reverse_string(path)
        
        # If the reverse path is not in final set, add the path.
        if revpath not in individual_paths:
            individual_paths.add(path)
    
    print 'individual', len(individual_paths)
    save_hamiltonian_paths(individual_paths, 'individual_paths.hpath')
    
    # After removing this degeneracy, there are 132 paths.
    # Each of these paths can start at any packaging signal position nominated.
    # For our analysis, we constrained the start and end of the paths to be at
    # the eight vertices nearest the maturation protein of MS2.
    
    # This meant the total realisation of paths on the sphere was 132*5*8=5,280
    
    # Note that this includes the fact that in our analysis, 5' and 3'
    # directions are not distinguishable, and also the trimming of start/end of
    # the path, discussed earlier.
    
### MAIN

if __name__ == '__main__':
    # Main primary forms a library from the generated protein positions
    
    main_primary()
    
    # Main secondary calculates the same library from lists of cycles and
    # pseudo cycles for MS2 (of length 60 & 59), discussed previously in the
    # literature. Included here for comparison purposes.
    
    #main_secondary()

### ENDS