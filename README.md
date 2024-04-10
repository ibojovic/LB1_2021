# **I Protein Structure Superimposition**
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

# **II Protein Sequence Alignment Analysis**

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
An example of the output: *cytc_aln.clw*

**parse_aln.py** script that analyzes the multiple sequence alignment and calculates for each position the most abundant residue with its frequency and the information entropy (S). Where S is

>$S=\sum_{i=0}^{20} -p_ilog(p_i)$

and *p* are the frequecnies of the amino acids.  The program has a function the parse the alignment file and a function that calculates the profile in each position of the alignment.

`!python3 parse_aln.py ../data/cytc_aln.clw`

The output of the program returns for each position the most frequent amino acid, the entropy, the frequency and the number of gaps.

Sorting this file you can detect the most conserved sites among the five proteins.

## **Results of a Blast search**

blast+ from the ftp website:

*   [Linux x64](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.11.0/ncbi-blast-2.11.0+-x64-linux.tar.gz)
*   [MacOSX](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.11.0/ncbi-blast-2.11.0+-x64-macosx.tar.gz)

Untar the file and download swissprot from the UniProt website

```bash
tar -xzvf ncbi-blast-2.11.0+-x64-[ARCH].tar.gz # Where [ARCH] is linux or macosx
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
```

Unzip the uniprot_sprot.fasta.gz file and format the database using makeblastdb
```bash
gunzip uniprot_sprot.fasta.gz
ncbi-blast-2.11.0+/bin/makeblastdb -in uniprot_sprot.fasta -input_type fasta -dbtype prot
```
Run the blast search on the Cytochrome C from Arabidopsis (Q93VA3) selecting sequences with E-value < 0.001.

`ncbi-blast-2.11.0+/bin/blastp -query Q93VA3.fasta -db uniprot_sprot.fasta -evalue 0.001 -out Q93VA3.blast -outfmt 7`

Parsing the output file to get the list of identifieer of the possible homolog proteins
```bash
# Fileter the blast output checking the e-value threshold ($11<=1e-2)
grep -v "^#" Q93VA3.blast |awk '{split($2,a,"|"); if ($11<=1e-3) print a[2]}' >Q93VA3.homolog

# This command do not filter the output of blast
grep -v "^#" Q93VA3.blast |cut -f 2 |cut -d "|" -f 2
```
Download the sequences of the homolog excluding the first element in the list (the query sequence (Q93VA3) previosly downloaded) Build a unique file containing all the sequence inluding Q93VA3
```bash
for i in `cat Q93VA3.homolog |tail -n +2`
do
   wget https://www.uniprot.org/uniprot/$i.fasta
done

for i in `cat Q93VA3.homolog`
do
   cat $i.fasta
done > cyt5_swiss.fasta
```

```bash
muscle[VERSION] -in cyt5_swiss.fasta -out cyt5_swiss.aln -clw 
```

```bash
!python3 parse_aln.py ../data/cyt5_swiss.aln P99999 |awk '{if ($3==1 && $5==0) print $0}'
```


# **III Probabilistic Models**

Transition probabilities that maximize the probability of the sequence of events

The script from the impot sequence generates the alphabet with the function get_alphabet.

The get_tmatrix funtion starting from position 0 and ending at the position n-1 counts all the transitions bewteen state i and i+1. In the final step of the function the matrix is normalized by line.

The normalization corrseponds to set to 1 the sum of probabilities of trantion from one state.

The function get_shuffle is written to generate random shuffled sequence and show that their probabilities are in general lower than the probability of the input sequence.

```bash
!python3 markov_model.py CCCFFCRRRCCSSSSFSFRRFFSSFSCSRRS
```

## **CpG Island HMM**
![image](https://github.com/ibojovic/protein_structure/assets/62520977/4395344f-664e-41dc-9672-bec97f18bc0d)

 Training sequence used here is the human chromosome 21 downloaded from the UCSC genome browser( **chr21.fa.gz** file). For the GpC Island annotation refer to the **cpgIslandExt.txt.gz** file.

 The script **cpg-hmm.py** takes in input bedfile extracted from the cpgIslandExt.txt and the sequence of the chromosome 21 and calculates the transition and emission probabilities for each state.

 The first 3 functions of the scripr parse the bedfile and sequfile and generated 2 strings containing the DNA sequence and the CpG state. Only standard nucleic acids are considered for the statistics.

This calculation is performed on the first 50 lines of the cpgIslandExt.txt file and the first 10 million nucleotides of chromosome 21. The pickle library is used to save the resulting HMM transission and probabilities on a file ('../data/cpg-hmm.pik') for future calculations.

```bash
python cpg-hmm.py ../data/chr21_short.cpg ../data/chr21_short.fa
```

**forward-hmm.py** script calculates the probability of the sequence giving the model using the forward algorithm.

```bash
python forward-hmm.py ACGTG ../data/cpg-hmm.pik
```










