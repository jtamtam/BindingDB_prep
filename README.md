# BindingDB_Preprocessing
Flowchart of the tasks:

<img width="526" alt="image" src="https://user-images.githubusercontent.com/95087658/171919481-de060115-901a-4938-8f92-56dba0517087.png">

## Advice
To run the different files I suggest the creation of a local environment and you will need to install ncbi-blast-2.13.0+.
The blastp exec and the makeblastdb exec must be in the /bin of the local environment. 
For the different programms to run correctly and without erros, they need to be launched in a specific order:
1. Itial_format.py
2. BLAST.py
3. Read_BLAST.py
4. Tanimoto.py
5. Structures.py

## Initial_format.py
This file takes BindingDB dataset, split it in subsets. There are 18 subsets based on (ki,kd,ic50)(train,test,valid)(with_duplicates, without_duplicates), which are saved in the "BindingDB_Preprocessing/data" folder
It produces FASTA files. The 9 files(without_duplicates) from the BindingDB dataset are saved in "BindingDB_Preprocessing/Fasta_files/bindingDB". It also produces a FASTA file from the general PDB, which is saved in "BindingDB_Preprocessing/Fasta_files/pdb".

## BLAST.py
This file executes a BLAST of the 9 files from the BindingDB (without_duplicates) against the pdb FASTA and vice-versa. The 18 resulting files are stored in "BindingDB_Preprocessing/Fasta_files/BLAST_Output"
(takes approx. 15 hours)

## Read_BLAST.py
This file processes the results of the BLAST and produces 9 dataframes with the following columns:
(Index) | Drug_ID | Drug | Target_ID | Target | Y(ic50/Kd/Ki) | Coverage | protein_id | protein_id_short |
------- | ------- | ---- | --------- | ------ | ------------- | -------- | ---------- | ---------------- |

The identity BLAST between the target and the protein is necessarily 100%

## Tanimoto.py
This file processes the results from Read_BLAST.py and compares the drugs SMILE with the ligands of the proteins in the PDB. It uses the Tanimoto calculation to assess a similarity coefficient.
It produces 9 final dataframes:
* train_ic50_final.csv
* test_ic50_final.csv
* valid_ic50_final.csv
* train_kd_final.csv
* test_kd_final.csv
* valid_kd_final.csv
* train_ki_final.csv
* test_ki_final.csv
* valid_ki_final.csv
         
Each of this dataframe has the following columns:
(Index)|Drug_ID|Drug|Target_ID|Target|Y(ic50/Kd/Ki)|Coverage|protein_id|protein_id_short|Ligand_Similarity|Ligand_ID|Ligand_Formula|
-------|-------|----|---------|------|-------------|--------|----------|----------------|-----------------|---------|--------------|

(takes approx. 3 hours)

## Structures.py
This files downloads the relevant proteins based  on the conditions of the flowchart. Here is the folder hierarchy in the Final_Data_Structures folder created in this programm.
Each folder contains the structures of the relevant proteins:
* train_kd_structures
  * coverage of 1
  * coverage > 0.8 and ligand similarity > 0.7
  * coverage > 0.8
  * coverage < 0.8 and ligand similarity  > 0.7
* test_kd_structures
  * coverage of 1
  * coverage > 0.8 and ligand similarity > 0.7
  * coverage > 0.8
  * coverage < 0.8 and ligand similarity  > 0.7
* valid_kd_structures
  * coverage of 1
  * coverage > 0.8 and ligand similarity > 0.7
  * coverage > 0.8
  * coverage < 0.8 and ligand similarity  > 0.7
* train_ic50_structures
  * coverage of 1
  * coverage > 0.8 and ligand similarity > 0.7
  * coverage > 0.8
  * coverage < 0.8 and ligand similarity  > 0.7
* test_ic50_structures
  * coverage of 1
  * coverage > 0.8 and ligand similarity > 0.7
  * coverage > 0.8
  * coverage < 0.8 and ligand similarity  > 0.7
* valid_ic50_structures
  * coverage of 1
  * coverage > 0.8 and ligand similarity > 0.7
  * coverage > 0.8
  * coverage < 0.8 and ligand similarity  > 0.7
* train_ki_structures
  * coverage of 1
  * coverage > 0.8 and ligand similarity > 0.7
  * coverage > 0.8
  * coverage < 0.8 and ligand similarity  > 0.7
* test_ki_structures
  * coverage of 1
  * coverage > 0.8 and ligand similarity > 0.7
  * coverage > 0.8
  * coverage < 0.8 and ligand similarity  > 0.7
* valid_ki_structures
  * coverage of 1
  * coverage > 0.8 and ligand similarity > 0.7
  * coverage > 0.8
  * coverage < 0.8 and ligand similarity  > 0.7

