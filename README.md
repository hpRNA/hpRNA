hpRNA README
============

### v. 15/01/2015
### James Geraets, University of York
### jg923@york.ac.uk

Code used to generate and interrogate Hamiltonian paths corresponding to bacteriophage MS2. For publication 2014. GLP v3 license applies: see separate file for information.

CONTENTS OF THIS FILE
---------------------
 
  * Introduction
  * Requirements
  * Usage
  * Examples
  * Configuration
  * License
  * Maintainance

INTRODUCTION
------------

This code, as supplied, generates Hamiltonian paths as models for the RNA organization in proximity to capsid, for the bacteriophage MS2. A script is also supplied that allows the user to apply constraints to the library of Hamiltonian paths, thus enabling the addition of constraints arising from further structural insights, such as restriction to circular Hamiltonian paths, that narrow down the possible path solutions. The software is provided "as is", but the author is willing to correspond to help anyone with interest to amend/adapt the code, to interrogate asymmetric structures of other viruses and determine possible RNA conformations. A link to a scientific manuscript describing the process will be added here upon publication. See the NOTES.md file for additional guidance, or the website at hpRNA.github.io.

REQUIREMENTS
------------

The software is coded in python 2.7, a scripting language that is very easy to use and adapt. More information on python can be found at python.org or greenteapress.com/thinkpython/ 

Required python modules are:

  * pycairo
  * matplotlib

Other required software:
  * cairo    cairographics.org

USAGE
-----

Please see NOTES.md file or website (hpRNA.github.io) for contextual guidance.

  * hpRNA_generate.py
  
Generate connected paths on a polyhedral cage

usage: 

    hpRNA_generate.py -c CONNECTIVITY -s START [-h] [-e END] [-r REQUIRE] [-p PRECLUDE] [-l LENGTH] [-i ITERATION] [-d DEGENERACY] [-b] [-o OUTPUT]

required arguments:

  -c CONNECTIVITY, --connectivity CONNECTIVITY
                        File CONNECTIVITY. Provide a neighbor connectivity
                        map. First column position linked to positions in
                        other columns. Number of links does not have to be
                        uniform.

  -s START, --start START
                        File START. Provide multiple starting positions (or
                        partial paths of same length), to begin each generated
                        path.

optional arguments:

  -h, --help            Provides these usage instructions, then exits.

  -e END, --end END     File END. Enforce the endings to the generated paths
                        that are considered complete. Is used to ensure
                        circularization.

  -r REQUIRE, --require REQUIRE
                        File REQUIRE. Position availability requiring
                        previously visited positions. A move to positions in
                        first column requires prior visitation to those in
                        subsequent columns.

  -p PRECLUDE, --preclude PRECLUDE
                        File PRECLUDE. Provide exclusion based on previously
                        visited positions. A move to positions in first column
                        cannot occur if those in subsequent columns have
                        previously been visited.

  -l LENGTH, --length LENGTH
                        File LENGTH. Length of paths to consider, otherwise
                        paths visiting every position in CONNECTIVITY are
                        assumed: i.e. Hamiltonian path.

  -i ITERATION, --iteration ITERATION
                        Int ITERATION. Resume previously started generation at
                        supplied ITERATION (corresponding to path length).

  -d DEGENERACY, --degeneracy DEGENERACY
                        File DEGENERACY. Rotational symmetry of polyhedron:
                        all moves. For example, if the cage has icosahedral
                        symmetry, there will be 12*5=60 identical symmetric
                        views. Each line of the file represents a rotation to
                        an identical view: in every row, the elements have the
                        same relationship to each other, but the rows begin at
                        different positions.

  -b, --both            Option. Paths are calculated both 5'-3' and 3'-5'.
                        This only will make a difference if --require or
                        --preclude are used. Requires --degeneracy.

  -o OUTPUT, --output OUTPUT
                        Directory OUTPUT. Choose output directory. Default
                        'paths'.

  * hpRNA_constrain.py
  
Realize and constrain connected paths mapping to a polyhedral cage.

usage:

    hpRNA_constrain.py -p PATHS [-h] [-x CONSTRAINTS] [-d DEGENERACY] [-r REALIZE] [-b] [-o OUTPUT] [-m] [-c CONNECTIVITY] [--ms2]

required arguments:

  -p PATHS, --paths PATHS
                        File PATHS. Provide paths to realize or constrain.

optional arguments:

  -h, --help            Provides these usage instructions, then exits.

  -x CONSTRAINTS, --constraints CONSTRAINTS
                        File CONSTRAINTS. Provide constraints for paths, i.e.
                        edges of the polyhedral cage that are either present
                        (1) or not present (0).

  -d DEGENERACY, --degeneracy DEGENERACY
                        File DEGENERACY. Rotational symmetry of polyhedron:
                        all moves. For example, if the cage has icosahedral
                        symmetry, there will be 12*5=60 identical symmetric
                        views. Each line of the file represents a rotation to
                        an identical view: in every row, the elements have the
                        same relationship to each other, but the rows begin at
                        different positions.

  -r REALIZE, --realize REALIZE
                        File REALIZE. Points to realize the paths from.
                        Generated paths start from a small subset of points,
                        specified in the START file. Realize creates copies of
                        these general paths, to the symmetric frames of points
                        given in the REALIZE file, with respect to the frame
                        of the initial path being 'a'. Requires --degeneracy.

  -b, --backwards       Option. Realize in both directions: paths are reversed
                        additionally. Requires --realize.

  -o OUTPUT, --output OUTPUT
                        Directory OUTPUT. Choose output directory. Default
                        'paths'.

  -m, --moves           Option. Abstracts to a numbered move view, suitable
                        for simple symmetric cages. Numbered moves are
                        allocated from CONNECTIVITY file, thus correct
                        ordering of row elements in CONNECTIVITY file is
                        essential. Requires --connectivity.

  -c CONNECTIVITY, --connectivity CONNECTIVITY
                        File CONNECTIVITY. Provide a neighbor connectivity
                        map. First column position linked to positions in
                        other columns. Number of links does not have to be
                        uniform.

  --ms2                 Option. Additional analysis to compare to published
                        example of bacteriophage ms2. Provides graphical
                        output. Requires --connectivity.

EXAMPLES
--------

example_1:

    ./hpRNA_generate.py -c example_1/connectivity.txt -s example_1/start.txt -e example_1/end.txt -o example_1

This would generate Hamiltonian paths for the ms2 example. Neighbor map corresponds (example_1/connectivity.txt) shows how the positions are linked. The paths start at position 'a' and end at postions 'd','u','N','t' which are at the same five-fold vertex as 'a' (see geometry_guide.png). Paths saved to example_1 folder. This script can run *very* slowly. A sped up version is provided as example 4.

example_2:

    ./hpRNA_generate.py -c example_2/connectivity.txt -s example_2/start.txt -l example_2/length.txt -o example_2 -b -d example_2/degeneracy.txt -r example_2/require.txt

This would generate connected paths corresponding to neighbor map in connectivity.txt (ms2), starting with 'ac', and of lengths 12 and 8. These would not be Hamiltonian paths. The paths are extended in each direction (using the both flag - compare with and without). Additionally, the entire first vertex requires nucleating before move to another position (see require.txt, every position requires the positions of the first vertex filled). Paths saved to example_2 folder.

example_3:

    ./hpRNA_generate.py -c example_3/connectivity.txt -s example_3/start.txt -o example_3

Generate connected paths on a dodecahedron, starting with 'ab'. The connectivity map in this folder is different, corresponding to the dodecahedral geometry. Paths saved to example_3 folder.

example_4:

    ./hpRNA_generate.py -c example_4/connectivity.txt -s example_4/start.txt -e example_4/end.txt -o example_4

Sped up version of example_1. The paths start with strings 'adc' or 'dab' and end with 'uNt' or'tNu' which are at the same five-fold vertex as 'a' and 'd' (see geometry_guide.png). Note that this will not create all paths from example_1, but a subset. However, the following script (hpRNA_constrain.py) generalizes the paths for ms2 using the --ms2 flag. Paths saved to example_4 folder.

example_5:

    ./hpRNA_constrain.py -p example_5/paths_out.txt -r example_5/realize.txt -d example_5/degeneracy.txt -m -o example_5

Realize instances of the 132 Hamiltonian paths for ms2, on 8 of the 12 vertices (proteins specified in realize.txt file). Gives 5280 paths to test against constraints formulated from analysis of the tomographic data.

example_6:

    ./hpRNA_constrain.py -p example_6/paths_out_realized.txt --ms2 -x example_6/constrain.txt -c example_6/connectivity.txt -o example_6

Constraints deriving from tomographic data are applied to the 5280 paths realized in example_5. This results in 5 possible results. Constraints are connections between positions, and are marked in the constrain.txt file with (1) indicating must be occupied and (0) indicating must not be occupied. Remaining edges are free to be either occupied or unoccupied. If output is below 20 paths, then graphical representations (corresponding to geometry_guide.png) are drawn. Here, green and red dashed refer to constraints from constrain.txt, with green indicating occupied constraints and red dashed indicating non-occupied constraints. The inferred paths are given in black. The --ms2 tag has also cleaved the start and end of paths around the starting/ending vertex.

CONFIGURATION
-------------

Given the above guide, and comments in the code, the program can be amended/adjusted to calculate Hamiltonian or non-Hamiltonian connected paths for a wide variety of scenarios, and additional constraints can be included if desired by the user. For support/collaboration, please contact jg923@york.ac.uk

LICENSE
-------

License applies to software. Please see separate license file.

MAINTAINANCE
------------

This software is currently maintained by James Geraets, University of York. Support is envisaged for at least a couple of years. Comments and queries should be addressed by email or to:

James Geraets
RCH\321 YCCSA
University of York
York
YO10 5GE
United Kingdom
