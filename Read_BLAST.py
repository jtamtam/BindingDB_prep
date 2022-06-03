import pandas as pd
import numpy as np
from tqdm import tqdm
import os

train_ic50__no_dupl_blast=pd.read_csv("Fasta_files/BLAST_Output/pdb_vs_train_ic50_no_dupl_blast.out", sep="\t", header=None)
test_ic50_no_dupl_blast=pd.read_csv("Fasta_Files/BLAST_Output/pdb_vs_test_ic50_no_dupl_blast.out", sep="\t", header=None)
valid_ic50_no_dupl_blast=pd.read_csv("Fasta_Files/BLAST_Output/pdb_vs_valid_ic50_no_dupl_blast.out", sep="\t", header=None)

train_ki__no_dupl_blast=pd.read_csv("Fasta_Files/BLAST_Output/pdb_vs_train_ki_no_dupl_blast.out", sep="\t", header=None)
test_ki_no_dupl_blast=pd.read_csv("Fasta_Files/BLAST_Output/pdb_vs_test_ki_no_dupl_blast.out", sep="\t", header=None)
valid_ki_no_dupl_blast=pd.read_csv("Fasta_Files/BLAST_Output/pdb_vs_valid_ki_no_dupl_blast.out", sep="\t", header=None)

train_kd__no_dupl_blast=pd.read_csv("Fasta_Files/BLAST_Output/pdb_vs_train_kd_no_dupl_blast.out", sep="\t", header=None)
test_kd_no_dupl_blast=pd.read_csv("Fasta_Files/BLAST_Output/pdb_vs_test_kd_no_dupl_blast.out", sep="\t", header=None)
valid_kd_no_dupl_blast=pd.read_csv("Fasta_Files/BLAST_Output/pdb_vs_valid_kd_no_dupl_blast.out", sep="\t", header=None)

train_ic50_dupl=pd.read_csv("data/train_ic50_dupl.csv",header=0,index_col=0)
test_ic50_dupl=pd.read_csv("data/test_ic50_dupl.csv",header=0,index_col=0)
valid_ic50_dupl=pd.read_csv("data/valid_ic50_dupl.csv",header=0,index_col=0)

train_ki_dupl=pd.read_csv("data/train_ki_dupl.csv",header=0,index_col=0)
test_ki_dupl=pd.read_csv("data/test_ki_dupl.csv",header=0,index_col=0)
valid_ki_dupl=pd.read_csv("data/valid_ki_dupl.csv",header=0,index_col=0)

train_kd_dupl=pd.read_csv("data/train_kd_dupl.csv",header=0,index_col=0)
test_kd_dupl=pd.read_csv("data/test_kd_dupl.csv",header=0,index_col=0)
valid_kd_dupl=pd.read_csv("data/valid_kd_dupl.csv",header=0,index_col=0)

train_ic50=pd.DataFrame()
test_ic50=pd.DataFrame()
valid_ic50=pd.DataFrame()

train_ki=pd.DataFrame()
test_ki=pd.DataFrame()
valid_ki=pd.DataFrame()

train_kd=pd.DataFrame()
test_kd=pd.DataFrame()
valid_kd=pd.DataFrame()

array_of_array_of_data=([[train_ic50__no_dupl_blast,train_ic50_dupl,train_ic50],\
    [test_ic50_no_dupl_blast,test_ic50_dupl,test_ic50],\
    [valid_ic50_no_dupl_blast,valid_ic50_dupl,valid_ic50],\
    [train_ki__no_dupl_blast,train_ki_dupl,train_ki],\
    [test_ki_no_dupl_blast,test_ki_dupl,test_ki],\
    [valid_ki_no_dupl_blast,valid_ki_dupl,valid_ki],\
    [train_kd__no_dupl_blast,train_kd_dupl,train_kd],\
    [test_kd_no_dupl_blast,test_kd_dupl,test_kd],\
    [valid_kd_no_dupl_blast,valid_kd_dupl,valid_kd]])

#only keeping the rows where we have 100% identity
for aad in array_of_array_of_data:
    aad[0].columns=["query_id", "subject_id", "%identity", "alignment_len","mismatches","gap_opens","q_start","q_end","s_start","s_end","e_value","bit_score"]
    aad[0]=aad[0].loc[aad[0]['%identity'] >= 100-10e-9]

for aad in tqdm(array_of_array_of_data):
    coverage=None
    index_to_keep=[]
    coverages=[]
    for i,row in tqdm(aad[0].iterrows()):
        target=aad[1][aad[1].Target_ID==row.subject_id].Target.iloc[0]
        coverage=row.alignment_len/len(target)
        coverages.append(coverage)   
    aad[0]["coverage"]=coverages


for aad in tqdm(array_of_array_of_data):
    for i,d in aad[0].iterrows():
        df=aad[1][aad[1]["Target_ID"]==d.subject_id]
        df.insert(len(df.columns),"coverage",d.coverage)
        df.insert(len(df.columns),"protein_id",d.query_id)
        df.insert(len(df.columns),"protein_id_short",d.query_id[:4])
        aad[2]=pd.concat([aad[2],df])

if not os.path.exists("Before_Ligands"):
    os.makedirs("Before_Ligands")

array_of_array_of_data[0][2].to_csv("Before_Ligands/train_ic50.csv")
array_of_array_of_data[1][2].to_csv("Before_Ligands/test_ic50.csv")
array_of_array_of_data[2][2].to_csv("Before_Ligands/valid_ic50.csv")

array_of_array_of_data[3][2].to_csv("Before_Ligands/train_ki.csv")
array_of_array_of_data[4][2].to_csv("Before_Ligands/test_ki.csv")
array_of_array_of_data[5][2].to_csv("Before_Ligands/valid_ki.csv")

array_of_array_of_data[6][2].to_csv("Before_Ligands/train_kd.csv")
array_of_array_of_data[7][2].to_csv("Before_Ligands/test_kd.csv")
array_of_array_of_data[8][2].to_csv("Before_Ligands/valid_kd.csv")


