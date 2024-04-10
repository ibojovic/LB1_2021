# **Protein Structure Superimposition**
## **Script for superimposing two 3D structures.**
The script *super_pdb.py** should take as an input two PDB structures, the relative chains and their residue intervals. The script includes the following steps:

1.   Parsing the PDB file
2.   Selection of the subset of residues
3.   Calculation of the transformation that minimizes the RMSD

The last step is performed using the [biopython](https://biopython.org/) library. In particular, the module [SVDSuperimposer](https://biopython.org/DIST/docs/api/Bio.SVDSuperimposer-module.html) that allows to calculate the best transformation and the RMSD.

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

`!wget  https://files.rcsb.org/view/3O20.pdb -O data/3O20.pdb`

`!python script/super_pdb.py data/3O20.pdb data/3ZCF.pdb A A 10-60 10-60`

**WARNING**: This code is not safe in presence of Alternate locations in column 17 of the ATOM field.








