###**Protein Structure Superimposition**
##**Script for superimposing two 3D 
structures.**
The script  that should take in input two PDB structures, the relative chains and their residue intervals. The script has to include the following steps:

1.   Parse the PDB file
2.   Select the subset of residues
3.   Calculate the transformation that minimizues the RMSD

The last step is should be performed using the [biopython](https://biopython.org/) library. In particular we should use the module [SVDSuperimposer](https://biopython.org/DIST/docs/api/Bio.SVDSuperimposer-module.html) that allow to calculate the best transformation and the RMSD.


For addressing this task we build a basic python script with two functions. 

1.   Parse the PDB file, select the subset of residues
2.   Call SVDSuperimposer and return the transformation (Rotation matrix and Transaltion vector) and the RMSD.

As a reference for the selection of the coordinates in the PDB file consider the following fields in the [ATOM record](https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM):

*   ATOM: 1-4
*   ATOM TYPE: 13-16
*   RESIDUE TYPE: 18-20
*   CHAIN: 22
*   RESIDUE NUMBER: 23-26
*   X-COORD: 31-38
*   Y-COORD: 39-46
*   Z-COORD: 47-54

To download a PDB structure directly form the command line use the commad wget. 

!wget  https://files.rcsb.org/view/3O20.pdb -O data/3O20.pdb

!python script/super_pdb.py data/3O20.pdb data/3ZCF.pdb A A 10-60 10-60







