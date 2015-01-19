#!/usr/bin/env python
################################################################################
##                                                                            ##
##  hpRNA_constrain.py                                                SCRIPT  ##
##  ------------------------------------------------------------------------  ##
##                                                                            ##
##  v. 15/01/2015                                                             ##
##                                                                            ##
##  James Geraets, University of York                                         ##
##  jg923@york.ac.uk                                                          ##
##                                                                            ##
##  Tool for generating instances of Hamiltonian paths satisfying viral       ##
##  geometries, starting and ending at specified positions. These instances   ##
##  of paths can be subsequently constrained to only paths that satisfy       ##
##  particular cage edge utilization.                                         ##
##                                                                            ##
##  (C) University of York 2014                                               ##
##                                                                            ##
################################################################################

### MODULE IMPORTS

import os
import sys
import argparse
import numpy as np
import string

### FUNCTION DEFINITIONS

def translate(path, point, degeneracy):
    '''
    Translate the path to start at required starting point.
    
    '''
    trans = degeneracy[point]
    return ''.join(trans[a] for a in path)

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
    
def notation(path, connectivity):
    '''
    Change to move notation, using connectivity map to allocate moves.
    
    '''
    return ''.join(str(connectivity[a].index(b) + 1) for a, b in zip(path[:-1], path[1:]))
        
def comparison(input_paths, output_paths):
    '''
    Print a little comparison before and after constraints
    
    '''
    
    best_paths = count_best_ms2(output_paths)
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

def constrain(args):
    '''
    Function for constraining paths. A constraint file is loaded containing
    components of the paths, with specifiers to whether the components are
    present or not present. These paths are filtered against these constraints:
    only paths that pass all constraints are output. If ms2 is selected,
    special analysis will be run (including drawing of output paths).
    
    '''
    
    constrain_occ = []
    constrain_unocc = []
    
    for line in args.constraints:
        (edgeA, edgeB), boolean = line.strip().split()
        if int(boolean):
            constrain_occ.append((edgeA+edgeB, edgeB+edgeA))
        else:
            constrain_unocc.append((edgeA+edgeB, edgeB+edgeA))
    args.constraints.close()

    # Load in file with constrained edges
    # Sample constraint file provided as constraint.txt

    # Sift through paths, removing those that do not meet constraints
    hpath_input_name, hpath_input_extension = os.path.splitext(os.path.basename(args.paths.name))
    
    infile = args.paths
    
    if not args.ms2:
        outfile_n = os.path.join(args.output, hpath_input_name + '_constrained' + hpath_input_extension)
        outfile = open(outfile_n, 'w')
    
    if args.moves or args.ms2:
        connfile = args.connectivity
        connectivity = {}
        for line in connfile:
            line = line.strip().split()
            connectivity[line[0]] = line[1:]
        connfile.close()
        if not args.ms2:
            m_outfile_n = os.path.join(args.output, hpath_input_name + '_constrained_moves' + hpath_input_extension)
            m_outfile = open(m_outfile_n, 'w')
    
    ms2_input_paths = []
    ms2_output_paths = []
    incount = 0
    outcount = 0
    
    for hampath in infile:
        incount += 1
        hampath = hampath.strip()
        if args.ms2:
            ms2_input_paths.append((notation(hampath, connectivity), hampath))
        test = True
        for opt1, opt2 in constrain_occ:
            if (opt1 not in hampath) and (opt2 not in hampath):
                test = False
                break
        for opt1, opt2 in constrain_unocc:
            if (opt1 in hampath) or (opt2 in hampath):
                test = False
                break
        if test:
            outcount += 1
            if args.ms2:
                ms2_output_paths.append((notation(hampath, connectivity), hampath))
            elif args.moves:
                outfile.write(hampath + '\n')
                m_outfile.write(notation(hampath, connectivity) + '\n')
            else:
                outfile.write(hampath + '\n')
    if args.ms2:
        input_paths = []
        output_paths = []
        for n, h in ms2_input_paths:
            input_paths.append((n[n.index('1'):n.rindex('1') + 1], h[n.index('1'):n.rindex('1') + 2]))
        for n, h in ms2_output_paths:
            output_paths.append((n[n.index('1'):n.rindex('1') + 1], h[n.index('1'):n.rindex('1') + 2]))

        input_paths = list(set(input_paths))
        output_paths = upshift_ms2(list(set(output_paths)))
        
        if len(output_paths) < 20 and len(output_paths) > 0:
            # If few result paths, display and draw.
            display_solution_paths(input_paths, output_paths[:5], args, constrain_occ, constrain_unocc)
        elif len(output_paths) > 0:
            # Otherwise, just list them.
            comparison(input_paths, output_paths)
        else:
            # If no result paths, the constraints you have imposed mean none are feasible.
            print '\nNO SOLUTIONS'
               
    elif args.moves:
        infile.close()
        outfile.close()
        m_outfile.close()
    else:
        infile.close()
        outfile.close()
        
def display_solution_paths(input_paths, output_paths, args, constrain_occ, constrain_unocc):
    
    comparison(input_paths, output_paths)
    
    import hpRNA_ms2_draw
    
    hpath_input_name, hpath_input_extension = os.path.splitext(os.path.basename(args.paths.name))
    
    print 'SOLUTION PATHS\n'
    for hampath, proteins in output_paths:
        draw = [a+b for a, b in zip(proteins[:-1], proteins[1:]) if ((a+b not in [e for tupl in constrain_occ for e in tupl]) and (a+b not in [e for tupl in constrain_unocc for e in tupl]))]
        pngname = os.path.join(args.output, hpath_input_name + '_output_' + proteins + '.png')
        hpRNA_ms2_draw.hami_draw(constrain_occ, constrain_unocc, draw, pngname, hampath, proteins)
    
def count_best_ms2(output_paths):
    '''
    Count the remaining paths that correspond to the best Hamiltonian path.
    Either of the four orientations of it left by hamiltonian_path_generation.py
    
    '''
    best_out = []
    
    for best_p in permute_path('13133331333313312222133331222213333122122221312133331331'):
        best_out += filter(lambda x: x[0] == best_p, output_paths)

    return len(best_out)
    
def upshift_ms2(paths_in, instances=None):
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

def realize(args):
    '''
    Function for creating instances of generalized paths, starting at different
    points on the polyhedron.
    
    '''
    degenmatrix = np.loadtxt(args.degeneracy, dtype=str)

    startp = degenmatrix[0,0]
    translation = {}
    for i in range(degenmatrix.shape[0]):
        translation[degenmatrix[i,0]] = dict(zip(degenmatrix[0,:], degenmatrix[i,:]))
    
    realizefile = args.realize
    args.realize = []
    for line in realizefile:
        args.realize.append(line.strip())
    realizefile.close()

    hpath_input_name, hpath_input_extension = os.path.splitext(os.path.basename(args.paths.name))

    outfile_n = os.path.join(args.output, '.' + hpath_input_name + '_realized_tmp' + hpath_input_extension)
    prunedfile_n = os.path.join(args.output, hpath_input_name + '_realized' + hpath_input_extension)
    
    infile = args.paths
    outfile = open(outfile_n, 'w')
    
    if args.moves:
        connfile = args.connectivity
        connectivity = {}
        for line in connfile:
            line = line.strip().split()
            connectivity[line[0]] = line[1:]
        connfile.close()
        m_outfile_n = os.path.join(args.output, '.' + hpath_input_name + '_moves_realized_tmp' + hpath_input_extension)
        m_prunedfile_n = os.path.join(args.output, hpath_input_name + '_moves_realized' + hpath_input_extension)
        
        m_outfile = open(m_outfile_n, 'w')
        
    for hampath in infile:
        hampath = hampath.strip()
        for realize_point in args.realize:
            trans = translate(hampath, realize_point, translation)
            backwards = backwards_string(trans)
            outfile.write(trans + '\n')
            if args.backwards:
                outfile.write(backwards + '\n')
            if args.moves:
                m_outfile.write(notation(trans, connectivity) + '\n')
                if args.backwards:
                    m_outfile.write(notation(trans, connectivity) + '\n')
    
    infile.close()
    outfile.close()
    
    if args.moves:
        m_outfile.close()
    
    awk = "awk \'!seen[$0]++\' " + outfile_n + " > " + prunedfile_n
    os.system(awk)
    
    if args.moves:
        awk_m = "awk \'!seen[$0]++\' " + m_outfile_n + " > " + m_prunedfile_n
        os.system(awk_m)
            

### MAIN

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Realize and constrain connected paths mapping to a polyhedral cage.")
    parser.add_argument("-p", "--paths", help='File PATHS. Provide paths to realize or constrain.', type=file, required=True)
    parser.add_argument("-x", "--constraints", help='File CONSTRAINTS. Provide constraints for paths, i.e. edges of the polyhedral cage that are either present (1) or not present (0).', type=file)
    parser.add_argument("-d", "--degeneracy", help='File DEGENERACY. Rotational symmetry of polyhedron: all moves. For example, if the cage has icosahedral symmetry, there will be 12*5=60 identical symmetric views. Each line of the file represents a rotation to an identical view: in every row, the elements have the same relationship to each other, but the rows begin at different positions.', type=file)
    parser.add_argument("-r", "--realize", help='File REALIZE. Points to realize the paths from. Generated paths start from a small subset of points, specified in the START file. Realize creates copies of these general paths, to the symmetric frames of points given in the REALIZE file, with respect to the frame of the initial path being \'a\'. Requires --degeneracy.', type=file)
    parser.add_argument("-b", "--backwards", help='Option. Realize in both directions: paths are reversed additionally. Requires --realize.', action='store_true')
    parser.add_argument("-o", "--output", help='Directory OUTPUT. Choose output directory. Default \'paths\'.', default='paths')
    parser.add_argument("-m", "--moves", help='Option. Abstracts to a numbered move view, suitable for simple symmetric cages. Numbered moves are allocated from CONNECTIVITY file, thus correct ordering of row elements in CONNECTIVITY file is essential. Requires --connectivity.', action='store_true')
    parser.add_argument("-c", "--connectivity", help='File CONNECTIVITY. Provide a neighbor connectivity map. First column position linked to positions in other columns. Number of links does not have to be uniform.', type=file)
    parser.add_argument("--ms2", help='Option. Additional analysis of constraints to compare to published example of bacteriophage ms2. Provides graphical output. Requires --connectivity and --constraints.', action='store_true')
    args = parser.parse_args()
    
    if (not args.output) and args.realize:
        args.output = 'realized'
    elif (not args.output) and (not args.realize):
        args.output = 'constrained'

    if not os.path.exists(args.output):
        try:
            os.makedirs(args.output)
        except:
            parser.error("--output directory error.")
            
    if args.realize and args.degeneracy is None:
        parser.error("--realize requires --degeneracy.")
    if args.backwards and args.realize is None:
        parser.error("--backwards requires --realize.")
    if args.moves and args.connectivity is None:
        parser.error("--moves requires --connectivity.")
    if args.ms2 and (args.connectivity is None or args.constraints is None):
        parser.error("--ms2 requires --connectivity and --constraints.")

    if args.realize:
        realize(args)
    elif args.constraints:
        constrain(args)
    else:
        parser.error("either --realize or --constraints is required")

### ENDS