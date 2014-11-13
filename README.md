hpRNA README
============

### v. 08/10/2014
### James Geraets, University of York
### jg923@york.ac.uk

Code used to generate and interrogate Hamiltonian paths corresponding to bacteriophage MS2. For publication 2014. GLP v3 license applies: see separate file for information.



CONTENTS OF THIS FILE

---------------------

  
 
  * Introduction

  * Requirements

  * Usage

  * Configuration

  * License

  * Maintainance




INTRODUCTION

------------



This code, as supplied, generates Hamiltonian paths that correspond to the 
possible RNA conformations in proximity to capsid, for the bacteriophage MS2.
A script is also supplied to apply constraints to the Hamiltonian paths due to
structural considerations, narrowing down the possible path solutions. The
software is provided "as is", but the author is willing to correspond to help
anyone with interest to amend/adapt the code, to interrogate asymmetric 
structures of other viruses and determine possible RNA conformations.




REQUIREMENTS

------------



The software is coded in python 2.7, a scripting language that is very easy to
use and adapt. More information on python can be found at python.org or 
greenteapress.com/thinkpython/

Required python modules are:


  * pycairo

  * matplotlib



Other required software:


  * cairo    cairographics.org




USAGE

-----



Step-by-step guide to recreating results for MS2.



### Generating Hamiltonian paths



Using hamiltonian_path_generate.py you can generate Hamiltonian paths that 
relate to the geometry specified within the script. The geometric cage is a 
result of ssRNA connecting all packaging signal (PS) positions on the 60
 heterodimers in the MS2 capsid. Neighbour maps, which contain detail about 
which PS positions are able to be moved between when realising the geometrical
 constraints, are depended on by the script.

These are set up for MS2 (or a
related phage such as GA), or another of the Leviviridae. 

However, the neighbour maps and general geometry in the script could be 
substituted, so to represent another virus with a different geometry. Other 
alterations would have to be made in following steps in the code. For example, 
for MS2, paths of length 59 are implicitly required to traverse the entire 
capsid. However, there are restrictions about the start/end vertex (which is 
designated as the same due to 5' and 3' ends binding to MP). The path cannot 
visit this vertex except at start/end (unlike other vertices). Other viruses 
will have different special cases, which requires careful consideration of the
literature.



We have utilised several steps to speed up the derivation of the Hamiltonian 
paths. However, for clarity, many of these have been rolled back for this open-
source code. Yet many special cases for MS2 have been coded, and must be 
considered when adapting the code: for example the binding of 5' and 3' ends of 
RNA to MP, which circularizes the RNA, and restricts the start/end vertex to near 
the MP.



Paths are calculated on a virus, starting at a certain protein, `a`, from the
list of proteins, which are layed out as shown in "geometry_guide.png":
 `abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSUVTWXYZ01234567`

Due to the restrictions on the primary vertex, we can reduce the search space by 
symmetry considerations. We can ignore "short" (C5 move) edges of the 
polyhedron: short moves around the five-fold primary vertex. This reduces the 
path to 57 moves, the first and last of which are both necessarily "long" (DS 
move) moves. 

Additionally, we reduce the search space by fixing a direction of the second 
move. This cannot be "long" (DS) as this would visit a previously visited heterodimer. Thus we arbitrarily pick a single direction. A mirror image path can be
created afterwards to recreate the ignored direction. 

Once paths have been generated on the protein neighbour map, they can be 
generalised into moves, i.e. described as the order of mapping operations 
between proteins, rather than the proteins themselves. This allows an easier way 
to recreate the paths starting from any given protein. Also, calculation of 
symmetric and mirror paths are much easier when the protein positions are 
disregarded: these degenerate paths can be calculated with simple string 
manipulations. This is all undertaken in the script 
hamiltonian_paths_multiply.py




### Realising Hamiltonian paths over the virus



Given the calculated generalised paths, we have to realise these paths starting 
and finishing at defined points on the viral capsid. These starting points
come from consideration of the binding of RNA 5' and 3' to the maturation
 protein, which localises the ends of the RNA in the vicinity of MP. We have 
called this "instantiating" individual paths, and this is undertaken in the 
script hamiltonian_path_instantiate.py 

Starting points are referred to by the same protein labelling map utilised in 
the first step:
 `abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSUVTWXYZ01234567`

### Difference maps


The analysis is based on a sub-tomographically-averaged, asymmetric structure of
 MS2 [Dent et al., Structure 2013]. The asymmetric structure was obtained by
 imaging mature MS2 bound to its natural receptor, the F-pilus of E. coli. The 
data can be found in the Electron Microscopy DataBank, with tag EMD-2365. The X-ray structure of the protein capsid of MS2 can be found at the Protein 
Data Bank, with PDBID 2MS2. This protein structure was filtered to 39AA 
resolution to match the EM data. Then the pixel size and orientation of the two 
maps are made equivalent by trilinear interpolation of the reduced-resolution 
X-ray structure. A contour mask of 0.5sigma is used to sample the low-resolution
 map and used to eliminate the protein density. 
A similar difference map is created between an icosahedrally-averaged map 
[Toropova et al., J Mol Biol 2008] and the resampled filtered protein, yielding
 a symmetric cage of RNA with a polyhedral shape.

 For both maps, the outer shell of RNA in proximity to capsid is isolated by
icosahedral masking with vertex radii of 80AA and 120AA.

### Mapping data onto the geometric model

The polyhedral density was partitioned into segments attributed to the edges of
the polyhedral cage of RNA seen in the symmetric map. Pixels from the asymmetric 
RNA map can be associated with defined segments on the polyhedral shell, and 
each connection thus has a density profile associated with it. 

Some segments were not used further in the analysis if they had not sampled 
many pixels, or were in parts of the tomogram where features would obscure the 
RNA density. In particular, connections adjacent to the MP/pilus are discarded 
as they may contain unmasked MP density.

### The density profiles


The densities of the connections were compared, and statements about the 
likelihood of RNA occupancy were made about these connections. These were noted 
as constraints. 

The connections were named in the following way. We had labelled the vertices as
these are the easiest to recognise in three-dimensional imaging programs. We 
have included a mapping of the vertices to the proteins, by way of the PS 
positions being labelled as proteins with the following letters and numbers:
 `abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSUVTWXYZ01234567`


"geometry_guide.png" shows the relation between the heterodimers and vertices. 
The maturation protein has been (arbitrarily) chosen to map to the homodimer on
the two-fold axis between vertices `5` and `6`. This is marked on the path outputs.




### Constraining paths



Constraints were saved in a file. An example file is ms2.constraints. In the 
first column, and `L` or `S` represents whether the constraint set is a long or 
short edge. The next 0 or 1 sets the edge either to unoccupied or occupied, 
respectively. Finally, two columns represented the name of the edge
 corresponding to the protein naming.

 For example, in the constraints file, a row of `L 1 e f` means that the long edge between the heterodimers `e` and `f` is constrained to 
be occupied in the analysis. Similarly `S 0 T Q` 
means the short edge between the heterodimers `Q` and `T` is constrained to be 
non-occupied in the analysis. Note that in the provided analysis, we were
unable to make use of any short-edge data, as the resolution of the asymmetric 
tomogram was too low. However, this did not prevent a good result using only 
long-edge data.

Paths corresponding to a constraints file can be found by running the 
hamiltonian_path_constrain.py script:
 `python <pick file of input paths> <constraints file>`




### Drawing Hamiltonian path solutions



A library for outputting figures of the solutions in 2d is given as 
hamiltonian_path_draw.py




CONFIGURATION

-------------



Given the above guide, and comments in the code, the program can be amended/
adjusted to calculate Hamiltonian paths for a wide variety of scenarios, and
 apply restrictions on these paths for the user. For support/collaboration,
 please contact jg923@york.ac.uk




LICENSE

-------



License applies to software. Please see separate license file.




MAINTAINANCE

------------


This software is currently maintained by James Geraets, University of York. 
Support is envisaged for at least a couple of years. Comments and queries should 
be addressed by email or to:


James Geraets  

RCH\321 YCCSA  

University of York  

York
YO10 5GE  

United Kingdom

