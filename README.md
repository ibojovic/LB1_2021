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
# **IV BLAST-based Annotation**
## **BLAST-based method for protein annotation.**

BLAST-based method for the annotation of BPTI/Kinutz domain.
Collection of a positive set of proteins containing the BPTI/Kinutz domain from UniProt and a negative set from the same database.
The human and non-human proteins are devided to build a training and testing sets respectively.
**Selecting training set from UniProt:**

*   Positives: database:(type:pfam pf00014) organism:"Homo sapiens (Human) [9606]" AND reviewed:y

*   Negatives: NOT database:(type:pfam pf00014) AND reviewed:yes AND organism:"Homo sapiens (Human) [9606]"



**Selecting testing set from UniProt:**

*   Positives: database:(type:pfam pf00014) NOT organism:"Homo sapiens (Human) [9606] AND reviewed:y

*   Negatives: NOT database:(type:pfam pf00014) AND reviewed:yes NOT organism:"Homo sapiens (Human) [9606]"

The composition of the sets is the following:

*  Human BPTI/Kunitz:  18
*  BPTI/Kunitz Non Human:  341
*  Not BPTI/Kunitz Human: 20347

For the negatives set is divided in two sets (Training and Testing) selecting the first 10000 proteins as training and the remaining ones as Testing. 

```bash
head -n 109074 Human_NotPF00014.fasta  >Human_NotPF00014_Training.fasta
tail -n +109075 Human_NotPF00014.fasta  >Human_NotPF00014_Testing.fasta
```
Format out positive training set of human BPTI/Kinutz proteins:
```bash
formatdb -i Human_PF00014.fasta
```
```bash
grep ">" Human_NotPF00014.fasta |head -n 10000 |tail -n 1
```

Run BLAST on the Human_PF00014.fasta dataset as a consistency test. 
```bash
# For checking purposes
blastpgp -i Human_PF00014.fasta -d Human_PF00014.fasta -o Human_PF00014.bl8 -m 8

# sort based on the 11 column to check the highst e-value
sort -grk 11 Human_PF00014.bl8| head -n 1
```
```bash
# Run on negatives
blastpgp -i Human_NotPF00014_Training.fasta -d Human_PF00014.fasta -o Human_NotPF00014_Training.bl8 -m 8

# sort the output to check the sequence with lowest e-value
sort -gk 11 Human_NotPF00014_Training.bl8 |head -n 2 |awk '{print $1,$2,$11}'
```
```bash
# Run on positive testing set
blastpgp -i NotHuman_PF00014.fasta -d Human_PF00014.fasta -o NotHuman_PF00014.bl8 -m 8

# Run on negative testing set
blastpgp -i Human_NotPF00014_Testing.fasta -d Human_PF00014.fasta -o Human_NotPF00014_Testing.bl8 -m 8
```
Rank the output on the posotive set of non Human BPTI/Kunitz proteins.
```bash
for i in `awk '{print $1}' NotHuman_PF00014.bl8 |sort -u `
do 
  grep $i NotHuman_PF00014.bl8 |sort -gk 11 |head -n 1
done > NotHuman_PF00014.bl8.best
```

Select an e-value of 0.001 as a classification threshold for the classification of BPTI/kunitz proteins. Accorddin to this assumption we calculated the number of proteins with maximu e-value below that threshold.
```bash
# Test on positives
awk '{if ($11<0.001) print $1}'  NotHuman_PF00014.bl8.best  |sort -u |wc -l

awk '{if ($11<0.001) print $1}'  Human_NotPF00014_Training.bl8  |sort -u |wc

awk '{if ($11<0.001) print $1}'  Human_NotPF00014_Tresting.bl8  |sort -u |wc -l
 ```

**performance.py** python script analyzes a file derived from the BLAST output.

The program will take in input a file with three columns representing the following data: 
* Protein ID
* The lowest e-value associated to each protein
* The class (0: NotBPTI/Kunitz 1: BPTI/Kunitz)

To run this progra we need the input file form the the BLAST output extracting the protein ID from the 1st column and the e-value from the 11th column and adding 0 or 1 depending on the protein files we are considering.
```bash
# Command for generating the positive set
awk '{split($1,v,'|'); print v[2],$11,1}' NotHuman_PF00014.bl8 > positive.txt

# Command for generating the negative Training set
awk '{split($1,v,'|'); print v[2],$11,0}' Human_NotPF00014_Training.bl8 > negative_train.txt

# Command for generating the negative Testing set
awk '{split($1,v,'|'); print v[2],$11,0}' Human_NotPF00014_Testing.bl8 > negative_test.txt
```
Merge positive and negatives in a unique file and give the resulting file in input to the python script above.
```bash
cat positive.txt negative_train.txt > train_data.txt
python performance.py train_data.txt

cat positive.txt negative_test.txt > test_data.txt
python performance.py test_data.txt
```

**WARNING**: BLAST do not return any output for sequence that match with an e-value higher tha 10 for this reason the total number of seuences in the pervious file will not sum up to the total number of proteins in the set.

Use the *comm* commad to determine the missing sequences and run again the script over the corrected files.

# **V Project Workflow**
## **Developemnet of a method for the detection BPTI/Kunitz domain in proteins.**

The workflow of the project is summarized as follows

1.  Selection of a representative set of protein structures from PDB.
2.  Multiple structural alignment with web available tools.
3.  Generation of the Hidden Markow Model for modeling BPTI/Kunitz domain.
4.  Selection of training and testing set from UniProt.
5.  Method optimization and assessment.

**Selection a representative set of protein structures from PDB.**

For this task we used the advanced search on the PDB web page. It is possible to use different contrains in the search. For the selection of the structures it is important to consider:

*   The resolution
*   The length of the protein
*   External identifiers for BPTI/Kunitz domain (e.g. PFAM, SCOP CATH)

When a list of PDB strutures is returned if it is feasible look to the proteins to check if they have a BPTI/Kunitz domain. 
From the PDB web page it is possible to download using the sequence reportthe protein sequence and check for the possible redundance in the set of proteins.
To check for the redundancy in the dataset you can use [*blastclast*](https://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/) algorithm. Blastclust takes as input  a fasta faile containing a set of sequences and cluster them according to the level of sequnece identity and coverage between the proteins. The option *-S* and *-L* for the program are used to control the percentage of sequence identity and the coverage respectively.

```bash
#The command sed and awk can be used to clean the csv file returned by the PDB.

# Replace inner comma
sed 's/, /_/g' filename 

# Remove spaces 
sed 's/ //g' filename

# split by comma
awk -F "," '{print $1,$2}'  filename

# After genereting a fasta file run blustclust with sequen identity of 95% and 90% coverage.
blastclust -i seqfile.fasta -o seqfile.clust -S 95 -L 0.9
```
Select a representative structure for each cluster considering the resolution of the structure and the length of the protein.

**Multiple structural alignment with web available tools.**

After collecting a clean set of protein structures calculate the multiple structure alignment using web available tools such as [PDBeFold](https://www.ebi.ac.uk/msd-srv/ssm/), [PROMALS3D](http://prodata.swmed.edu/promals3d) [mTM-align](https://yanglab.nankai.edu.cn/mTM-align/).
Check the returned alignment and look for the RMSD between pair of structures and the sequence identity. This information can be important to detect possible errors in the initial selection of the BPTI/Kunitz domain proteins.
Finally download the fasta file representing the multiple structure alignment.

**Generation of the Hidden Markow Model for modeling BPTI/Kunitz domain.**

When the alignemnt is returned it is importat to manually check the alignement to look for the conserved residues and the possible errors in the alignement. Ater that generate the HMM using *hmmbuild* command from [HMMER](http://hmmer.org/) package. 

```bash
# Generate the HMM
hmmbuild bpti-kunitz.hmm align-bpti-kunitz.fasta
```

**Selection of training and testing set from UniProt.**

This task is performed using the advanced search of UniProt database using as a key search the Pfam identifier of the BPTI/Kunitz domain. The search has to be restricted to the protein in SwissProt are manually annotated.
With this search  ~350 protein will be obtained containing BPTI/Kuntiz domain that will rapresente the positive set. All the remaining proteins in SwissProt can be used as negatives. 

For a fair test of  HMM model remove the sequences in the positive set that share high level of sequence identity with the protein structures collected for generating the model. This task can be performed using *blastpgp* with the option -m 8 checking for sequences with low e-value (column 11) and high sequence identity (column 3).

```bash
# Run trasforming the sequences of the selected structures in database
formatdb -i selected-bpti-kunits.fasta

# Run blastpgp to matche the positive set agains the selected bpti-kunitz.
blastpgp -i positives.fasta -d selected-bpti-kunits.fasta -m 8 -o positives.bl8
```

Using the output of blastpgp you can select the proteins that should be removed from the set of *positives.txt*

**Method optimization and assessment.**

When the final benchmark set has been collected a basic testing procedure consist in the implementation of a 2-fold cross-validation test. It consists in splitting positives and negatives in 2 subsets, optimizing the classification threshold on one subset and testing the performance on the other subset.

When the the dataset is splitted in two parts with the same size run the *hmmsearch* on the different subsets using following arguments: 

*   --max: turns off all the heuristics for cutting of distantly related proteins
*   --noali: exclude from the output the alignemnts 
*   --tblout: returns the output in tabular form
*   -Z: It is important for normalizing the e-value output 
*   --domZ: It is important for normalizing the domainn e-value output

An example of  command for macthing a set of sequences (clean_set.fasta) with the HMM model (bpti-kunitz.hmm) is the follwing:

```bash
# hmmsearch command matching the model (bpti-kunitz.hmm) with sequence in clean_set.fasta

hmmsearch --noali --max --tblout output.txt -Z 1 bpti-kunitz.hmm clean_set.fasta
```

The output file can be parsed considering the lines without the hash character (#) at the beginning selecting and the column 1 corresponding to the protein identifier and either column 5 or 8 that correspond to the fill-sequence and best-domain e-values respectively.

For the optimization and testing of the performace ckeck the hmmsearch output and include those proteins, presumably negatives, for which no match is returned. For the classification purposes assign that proteins and e-value greater or equal than the e-value theshold used for running hmmsearch.

When all the results for the subset are collected in files containg protein ID, e-value and protein class (0/1), we can use the **performace.py** script  to calculate the performance of our classification method at different e-value threshold.


























