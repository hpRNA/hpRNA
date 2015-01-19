#!/usr/bin/env python
################################################################################
##                                                                            ##
##  hpRNA_generate.py                                                 SCRIPT  ##
##  ------------------------------------------------------------------------  ##
##                                                                            ##
##  v. 15/01/2015                                                             ##
##                                                                            ##
##  James Geraets, University of York                                         ##
##  jg923@york.ac.uk                                                          ##
##                                                                            ##
##  Tool for generating a library of connected paths on a polyhedral cage,    ##
##  that can be used to interrogate asymmetric reconstructions of viruses.    ##
##                                                                            ##
##  (C) University of York 2014                                               ##
##                                                                            ##
################################################################################

### MODULE IMPORTS

import os
import sys
import subprocess
import numpy as np
import argparse
import shutil

### FUNCTION DEFINITIONS

def generate_paths(args):
    '''
    Main function to find Hamiltonian paths. Lengthens each input path by every
    single move possible, lengthening by one. Can be restarted if the run is
    incomplete, if given the last fully successful iteration as an argument.
    
    '''
        
    connfile = args.connectivity
    args.connectivity = {}
    for line in connfile:
        line = line.strip().split()
        args.connectivity[line[0]] = line[1:]
    connfile.close()
    
    firststart = args.start.readline().strip()
    args.start.seek(0)
    
    if not all(len(line.strip()) == len(firststart) for line in args.start):
        raise Exception("Not all start positions/paths are same length")
    args.start.seek(0)
        
    if args.iteration == None:
        args.iteration = len(firststart)
        shutil.copyfile(args.start.name, os.path.join(args.output, 'paths_%02i.txt' % (args.iteration,)))
        
    if args.degeneracy and args.both:
        degenfile = args.degeneracy
        degenmatrix = np.loadtxt(degenfile, dtype=str)
        degenfile.close()
        args.degeneracy = {}
        for i in range(degenmatrix.shape[0]):
            args.degeneracy[degenmatrix[i,0]] = dict(zip(degenmatrix[i,:], degenmatrix[0,:]))
    
    iteration = args.iteration
    
    if args.length:
        lengths = [int(a) for l in args.length for a in l.strip().split()]
        maxlength = max(lengths)
    else:
        lengths = [len(args.connectivity),]
        maxlength = len(args.connectivity)
    
    if args.require:
        req = dict((a.strip().split()[0], a.strip().split()[1:]) for a in args.require)
    else:
        req = None
    
    if args.preclude:
        pre = dict((a.strip().split()[0], a.strip().split()[1:]) for a in args.preclude)
    else:
        pre = None
    
    toolbar_width = maxlength - args.iteration
    
    # setup toolbar
    sys.stdout.write("[%s]" % (" " * toolbar_width))
    sys.stdout.flush()
    sys.stdout.write("\b" * (toolbar_width+1)) # return to start of line, after '['
        
    
    while iteration < maxlength:
        infile = open(os.path.join(args.output,'paths_%02i.txt' % (iteration,)), 'r')
        
        iteration += 1
        #print iteration
        
        if args.both == True:
            outfile = open(os.path.join(args.output,'.paths_%02i_tmp.txt' % (iteration,)), 'w')
            execute = execute_both
        else:
            outfile = open(os.path.join(args.output,'paths_%02i.txt' % (iteration,)), 'w')
            execute = execute_forward
        
        for i in infile:
            execute(i.strip(), args.connectivity, outfile, args.degeneracy, req, pre)
        
        infile.close()
        outfile.close()
        
        if args.both == True:
            awk = ("awk \'!seen[$0]++\' " +
                  os.path.join(args.output,'.paths_%02i_tmp.txt' % (iteration,)) + " > " +
                  os.path.join(args.output,'paths_%02i.txt' % (iteration,)))
            os.system(awk)
            
        sys.stdout.write("-")
        sys.stdout.flush()
    sys.stdout.write("\n")
    
    if args.end:
        endfile = args.end
        args.end = []
        for line in endfile:
            args.end.append(line.strip())
        endfile.close()
        
        outfile = open(os.path.join(args.output,'paths_out.txt'), 'w')
        
        for l in lengths:
            infile = open(os.path.join(args.output,'paths_%02i.txt' % (l,)), 'r')
        
            for i in infile:
                i=i.strip()
                for e in args.end:
                    if i[-len(e):] == e:
                        outfile.write(i+'\n')
            
            infile.close()
        outfile.close()
        
    else:
        outfile = open(os.path.join(args.output,'paths_out.txt'), 'w')
        
        for l in lengths:
            infile = open(os.path.join(args.output,'paths_%02i.txt' % (l,)), 'r')
        
            for i in infile:
                outfile.write(i)
            
            infile.close()
        outfile.close()


def execute_forward(i, connectivity, outfile, degeneracy, req, pre):
    '''
    Moves. Extend each path by one unit, attempting each possible move.
    
    '''
    
    for np in connectivity[i[-1]]:
        if np not in i:
            if not args.preclude or not any((p in i) for p in pre[np]):
                if not args.require or all((t in i) for t in req[np]):
                    outfile.write(i+np+'\n')


def execute_both(i, connectivity, outfile, degeneracy, req, pre):
    '''
    Moves from both 5' and 3' ends. Extend each path by one unit, attempting 
    each possible move.
    
    '''
    # Forwards
    for np in connectivity[i[-1]]:
        if np not in i:
            if not args.preclude or not any((p in i) for p in pre[np]):
                if not args.require or all((t in i) for t in req[np]):
                    outfile.write(i+np+'\n')
    # Backwards
    for np in connectivity[i[0]]:
        if np not in i:
            if not args.preclude or not any((p in i) for p in pre[np]):
                if not args.require or all((t in i) for t in req[np]):
                    outfile.write(rework(np+i, degeneracy)+'\n')

def rework(path, degeneracy):
    '''
    After backwards move, in 5' direction, the path can be redescribed without 
    loss of generality, to start at the same beginning point as all other 
    paths. This uses the degeneracy file supplied as a second optional argument 
    to the program.
    
    '''
    trans = degeneracy[path[0]]
    return ''.join(trans[a] for a in path)

### MAIN

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate connected paths on a polyhedral cage")
    parser.add_argument("-c", "--connectivity", help='File CONNECTIVITY. Provide a neighbor connectivity map. First column position linked to positions in other columns. Number of links does not have to be uniform.', type=file, required=True)
    parser.add_argument("-s", "--start", help='File START. Provide multiple starting positions (or partial paths of same length), to begin each generated path.', type=file, required=True)
    parser.add_argument("-e", "--end", help='File END. Enforce the endings to the generated paths that are considered complete. Is used to ensure circularization.', type=file)
    parser.add_argument("-r", "--require", help='File REQUIRE. Position availability requiring previously visited positions. A move to positions in first column requires prior visitation to those in subsequent columns.', type=file)
    parser.add_argument("-p", "--preclude", help='File PRECLUDE. Provide exclusion based on previously visited positions. A move to positions in first column cannot occur if those in subsequent columns have previously been visited.', type=file)
    parser.add_argument("-l", "--length", help='File LENGTH. Length of paths to consider, otherwise paths visiting every position in CONNECTIVITY are assumed: i.e. Hamiltonian path.', type=file)
    parser.add_argument("-i", "--iteration", help='Int ITERATION. Resume previously started generation at supplied ITERATION (corresponding to path length).', type=int)
    parser.add_argument("-d", "--degeneracy", help='File DEGENERACY. Rotational symmetry of polyhedron: all moves. For example, if the cage has icosahedral symmetry, there will be 12*5=60 identical symmetric views. Each line of the file represents a rotation to an identical view: in every row, the elements have the same relationship to each other, but the rows begin at different positions.', type=file)
    parser.add_argument("-b", "--both", help='Option. Paths are calculated both 5\'-3\' and 3\'-5\'. This only will make a difference if --require or --preclude are used. Requires --degeneracy.', action="store_true")
    parser.add_argument("-o", "--output", help='Directory OUTPUT. Choose output directory. Default \'paths\'.', default='paths')
    args = parser.parse_args()

    if not os.path.exists(args.output):
        try:
            os.makedirs(args.output)
        except:
            parser.error("--output directory error.")
    if args.both is True and args.degeneracy is None:
        parser.error("--both requires --degeneracy.")

    generate_paths(args)
    
## ENDS