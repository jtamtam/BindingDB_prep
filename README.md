# BindingDB_Preprocessing

## Advice
To run the different files I suggest the creation of a local environment and you will need to install ncbi-blast-2.13.0+.
The blastp exec and the makeblastdb exec must be in the /bin of the local environment. 

## Init_format.py
This file takes BindingDB dataset, split it in subsets. There are 18 subsets based on (ki,kd,ic50)(train,test,valid)(with_duplicates, without_duplicates), which are saved in the "BindingDB_Preprocessing/data" folder
It produces FASTA files. The 9 files(without_duplicates) from the BindingDB dataset are saved in "BindingDB_Preprocessing/Fasta_files/bindingDB". It also produces a FASTA file from the general PDB, which is saved in "BindingDB_Preprocessing/Fasta_files/pdb".

## BLAST.py
This file executes a BLAST of the 9 files from the BindingDB (without_duplicates) against the pdb FASTA and vice-versa. The 18 resulting files are stored in "BindingDB_Preprocessing/Fasta_files/BLAST_Output"

## Read_BLAST.py
This file processes the results of the BLAST and produces 9 dataframes with the following columns:
First Header  | Second Header
------------- | -------------

## Tanimoto.py

## Structures.py
