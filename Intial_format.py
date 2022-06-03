from tdc.multi_pred import DTI
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

import pandas as pd
import numpy as np
import os
import wget
import urllib.request
import gzip
import shutil

path="Fasta_files"
bindingDB_path=path+"/bindingDB"
pdb_path=path+"/PDB"
if not os.path.exists(path):
    os.makedirs(path)
if not os.path.exists(bindingDB_path):
    os.makedirs(bindingDB_path)
if not os.path.exists(pdb_path):
    os.makedirs(pdb_path)
if not os.path.exists("data"):
    os.makedirs(pdb_path)

data_kd = DTI(name = 'BindingDB_Kd')
data_ic50 = DTI(name = 'BindingDB_IC50')
data_ki = DTI(name = 'BindingDB_Ki')
split_kd = data_kd.get_split()
split_ic50 = data_ic50.get_split()
split_ki = data_ki.get_split()


train_kd = pd.DataFrame.from_dict(split_kd['train'])
train_kd_dupl=train_kd.dropna()
train_kd_dupl.to_csv("data/train_kd_dupl.csv")
train_kd_no_dupl=train_kd_dupl.drop_duplicates(subset=["Target_ID","Target"],inplace=False)
train_kd_no_dupl.to_csv("data/train_kd_no_dupl.csv")

valid_kd = pd.DataFrame.from_dict(split_kd['valid'])
valid_kd_dupl=valid_kd.dropna()
valid_kd_dupl.to_csv("data/valid_kd_dupl.csv")
valid_kd_no_dupl=valid_kd_dupl.drop_duplicates(subset=["Target_ID","Target"],inplace=False)
valid_kd_no_dupl.to_csv("data/valid_kd_no_dupl.csv")

test_kd = pd.DataFrame.from_dict(split_kd['test'])
test_kd_dupl=test_kd.dropna()
test_kd_dupl.to_csv("data/test_kd_dupl.csv")
test_kd_no_dupl=test_kd_dupl.drop_duplicates(subset=["Target_ID","Target"],inplace=False)
test_kd_no_dupl.to_csv("data/test_kd_no_dupl.csv")

train_ki = pd.DataFrame.from_dict(split_ki['train'])
train_ki_dupl=train_ki.dropna()
train_ki_dupl.to_csv("data/train_ki_dupl.csv")
train_ki_no_dupl=train_ki_dupl.drop_duplicates(subset=["Target_ID","Target"],inplace=False)
train_ki_no_dupl.to_csv("data/train_ki_no_dupl.csv")

valid_ki = pd.DataFrame.from_dict(split_ki['valid'])
valid_ki_dupl=valid_ki.dropna()
valid_ki_dupl.to_csv("data/valid_ki_dupl.csv")
valid_ki_no_dupl=valid_ki_dupl.drop_duplicates(subset=["Target_ID","Target"],inplace=False)
valid_ki_no_dupl.to_csv("data/valid_ki_no_dupl.csv")

test_ki = pd.DataFrame.from_dict(split_ki['test'])
test_ki_dupl=test_ki.dropna()
test_ki_dupl.to_csv("data/test_ki_dupl.csv")
test_ki_no_dupl=test_ki_dupl.drop_duplicates(subset=["Target_ID","Target"],inplace=False)
test_ki_no_dupl.to_csv("data/test_ki_no_dupl.csv")

train_ic50 = pd.DataFrame.from_dict(split_ic50['train'])
train_ic50_dupl=train_ic50.dropna()
train_ic50_dupl.to_csv("data/train_ic50_dupl.csv")
train_ic50_no_dupl=train_ic50_dupl.drop_duplicates(subset=["Target_ID","Target"],inplace=False)
train_ic50_no_dupl.to_csv("data/train_ic50_no_dupl.csv")

valid_ic50 = pd.DataFrame.from_dict(split_ic50['valid'])
valid_ic50_dupl=valid_ic50.dropna()
valid_ic50_dupl.to_csv("data/valid_ic50_dupl.csv")
valid_ic50_no_dupl=valid_ic50_dupl.drop_duplicates(subset=["Target_ID","Target"],inplace=False)
valid_ic50_no_dupl.to_csv("data/valid_ic50_no_dupl.csv")

test_ic50 = pd.DataFrame.from_dict(split_ic50['test'])
test_ic50_dupl=test_ic50.dropna()
test_ic50_dupl.to_csv("data/test_ic50_dupl.csv")
test_ic50_no_dupl=test_ic50_dupl.drop_duplicates(subset=["Target_ID","Target"],inplace=False)
test_ic50_no_dupl.to_csv("data/test_ic50_no_dupl.csv")




train_kd_no_dupl_fasta=[SeqRecord(Seq(targ),id=str(i)) for i,targ in zip(train_kd_no_dupl.Target_ID, train_kd_no_dupl.Target)]
SeqIO.write(train_kd_no_dupl_fasta, bindingDB_path+"/train_kd_no_dupl.fasta", "fasta")
valid_kd_no_dupl_fasta=[SeqRecord(Seq(targ),id=str(i)) for i,targ in zip(valid_kd_no_dupl.Target_ID, valid_kd_no_dupl.Target)]
SeqIO.write(valid_kd_no_dupl_fasta, bindingDB_path+"/valid_kd_no_dupl.fasta", "fasta")
test_kd_no_dupl_fasta=[SeqRecord(Seq(targ),id=str(i)) for i,targ in zip(test_kd_no_dupl.Target_ID, test_kd_no_dupl.Target)]
SeqIO.write(test_kd_no_dupl_fasta, bindingDB_path+"/test_kd_no_dupl.fasta", "fasta")

train_ki_no_dupl_fasta=[SeqRecord(Seq(targ),id=str(i)) for i,targ in zip(train_ki_no_dupl.Target_ID, train_ki_no_dupl.Target)]
SeqIO.write(train_ki_no_dupl_fasta, bindingDB_path+"/train_ki_no_dupl.fasta", "fasta")
valid_ki_no_dupl_fasta=[SeqRecord(Seq(targ),id=str(i)) for i,targ in zip(valid_ki_no_dupl.Target_ID, valid_ki_no_dupl.Target)]
SeqIO.write(valid_ki_no_dupl_fasta, bindingDB_path+"/valid_ki_no_dupl.fasta", "fasta")
test_ki_no_dupl_fasta=[SeqRecord(Seq(targ),id=str(i)) for i,targ in zip(test_ki_no_dupl.Target_ID, test_ki_no_dupl.Target)]
SeqIO.write(test_ki_no_dupl_fasta, bindingDB_path+"/test_ki_no_dupl.fasta", "fasta")

train_ic50_no_dupl_fasta=[SeqRecord(Seq(targ),id=str(i)) for i,targ in zip(train_ic50_no_dupl.Target_ID, train_ic50_no_dupl.Target)]
SeqIO.write(train_ic50_no_dupl_fasta, bindingDB_path+"/train_ic50_no_dupl.fasta", "fasta")
valid_ic50_no_dupl_fasta=[SeqRecord(Seq(targ),id=str(i)) for i,targ in zip(valid_ic50_no_dupl.Target_ID, valid_ic50_no_dupl.Target)]
SeqIO.write(valid_ic50_no_dupl_fasta, bindingDB_path+"/valid_ic50_no_dupl.fasta", "fasta")
test_ic50_no_dupl_fasta=[SeqRecord(Seq(targ),id=str(i)) for i,targ in zip(test_ic50_no_dupl.Target_ID, test_ic50_no_dupl.Target)]
SeqIO.write(test_ic50_no_dupl_fasta, bindingDB_path+"/test_ic50_no_dupl.fasta", "fasta")




if not os.path.exists("pdb_seqres.txt"):
    wget.download(url="https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz")
    with gzip.open('pdb_seqres.txt.gz', 'rb') as f_in:
        with open('pdb_seqres.txt', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

with open("pdb_seqres.txt") as f:
    pdb=f.read()
pdb=pdb.split("\n")


pdb_names=[]
pdb_id=[]
for i,data in enumerate(pdb[:-1]):
    if i%2==0:
        name=data[1:]
        pdb_names.append(name)
    else:
        pdb_id.append(data)
df_pdb=pd.DataFrame([pdb_names,pdb_id])
df_pdb=df_pdb.transpose()
df_pdb.columns=["names","target"]
df_pdb.drop_duplicates(subset="target",inplace=True)
#%%
pdb_file=[SeqRecord(Seq(frame.target),id=str(frame.names)) for _,frame in df_pdb.iterrows()]
SeqIO.write(pdb_file, pdb_path+"/pdb.fasta", "fasta")