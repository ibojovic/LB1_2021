# **Protein Structure Superimposition**
## **Script for superimposing two 3D structures.**
The script **super_pdb.py** should take as an input two PDB structures, the relative chains and their residue intervals. The script includes the following steps:

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

# **Protein Sequence Alignment Analysis**

## **Analysis of  a multuple sequence alignment to detect conserved sites.**

A unique fasta file containing all the sequences of the proteins downloaded from the UniProt Rest API 
```bash
for i in P99999 P00004 P0C0X8 P00091 Q93VA3
do
  wget https://www.uniprot.org/uniprot/$i.fasta
done
cat P99999.fasta P00004.fasta P0C0X8.fasta P00091.fasta Q93VA3.fasta > cytc_aln.clw
```

Calculate the multiple sequence alignment: [CLUSTALW](https://www.ebi.ac.uk/Tools/msa/clustalo/), [MUSCLE](https://www.ebi.ac.uk/Tools/msa/muscle/) and [TCOFFEE](https://www.ebi.ac.uk/Tools/msa/tcoffee/).

You should save the output in clustalw format that contains seperated columns representing the identifiers and the relative aligned sequences.

An example of the output: 

Now you should write a python script that analyzes the multiple sequence alignment and calculates for each position the most abundant residue with its frequencyy and the information entropy (S). Where S is

>$S=\sum_{i=0}^{20} -p_ilog(p_i)$

and *p* are the frequecnies of the amino acids.  The program will have a function the parse the alignment file and a function that calculates the profile in each position of the alignment.







