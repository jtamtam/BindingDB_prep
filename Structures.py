import pandas as pd
from Bio.PDB import PDBList
from os import path
import os
from tqdm import tqdm

if not os.path.exists("Final_Data_Structures"):
    os.makedirs("Final_Data_Structures/")

Path="Final_Data_Structures/"

train_kd_for_struct=pd.read_csv("Final_Data/train_kd_final.csv",index_col=0)
test_kd_for_struct=pd.read_csv("Final_Data/test_kd_final.csv",index_col=0)
valid_kd_for_struct=pd.read_csv("Final_Data/valid_kd_final.csv",index_col=0)

train_ic50_for_struct=pd.read_csv("Final_Data/train_ic50_final.csv",index_col=0)
test_ic50_for_struct=pd.read_csv("Final_Data/test_ic50_final.csv",index_col=0)
valid_ic50_for_struct=pd.read_csv("Final_Data/valid_ic50_final.csv",index_col=0)

train_ki_for_struct=pd.read_csv("Final_Data/train_ki_final.csv",index_col=0)
test_ki_for_struct=pd.read_csv("Final_Data/test_ki_final.csv",index_col=0)
valid_ki_for_struct=pd.read_csv("Final_Data/valid_ki_final.csv",index_col=0)

final_data=[
    [train_kd_for_struct,"train_kd_structures"],
    [test_kd_for_struct,"test_kd_structures"],
    [valid_kd_for_struct,"valid_kd_structures"],
    [train_ic50_for_struct,"train_ic50_structures"],
    [test_ic50_for_struct,"test_ic50_structures"],
    [valid_ic50_for_struct,"valid_ic50_structures"],
    [train_ki_for_struct,"train_ki_structures"],
    [test_ki_for_struct,"test_ki_structures"],
    [valid_ki_for_struct,"valid_ki_structures"],
]

pdbl = PDBList()

for fd in tqdm(final_data):

    if not os.path.exists(Path+fd[1]):
        os.makedirs(Path+fd[1])

    data_100_cov = fd[0].drop(fd[0][(fd[0]['coverage'] != 1)].index)
    list_100_cov=data_100_cov["protein_id_short"].tolist()
    data_100_cov.to_csv(Path+fd[1]+"/data_100_cov.csv")
    if not os.path.exists(Path+fd[1]+"/100_cov_struct"):
        os.makedirs(Path+fd[1]+"/100_cov_struct")
    pdbl.download_pdb_files(list_100_cov, pdir=Path+fd[1]+"/100_cov_struct", obsolete=False)

    data_80_cov_70_sim = fd[0].drop(fd[0][(fd[0]['coverage'] < 0.8)].index)
    data_80_cov_70_sim = data_80_cov_70_sim.drop(data_80_cov_70_sim[(data_80_cov_70_sim["Ligand_Similarity"] < 0.7)].index)
    data_80_cov_70_sim = data_80_cov_70_sim.drop(data_80_cov_70_sim[(data_80_cov_70_sim["coverage"] == 1)].index)
    list_80_cov_70_sim=data_80_cov_70_sim["protein_id_short"].tolist()
    data_80_cov_70_sim.to_csv(Path+fd[1]+"/data_80_cov_70_sim.csv")
    if not os.path.exists(Path+fd[1]+"/80_cov_70_sim_struct"):
        os.makedirs(Path+fd[1]+"/80_cov_70_sim_struct")
    pdbl.download_pdb_files(list_80_cov_70_sim, pdir=Path+fd[1]+"/80_cov_70_sim_struct", obsolete=False)

    data_80_cov = fd[0].drop(fd[0][(fd[0]['coverage'] < 0.8)].index)
    data_80_cov = data_80_cov.drop(data_80_cov[(data_80_cov['Ligand_Similarity'] > 0.7)].index)
    data_80_cov = data_80_cov.drop(data_80_cov[(data_80_cov['coverage'] == 1)].index)
    list_80_cov=data_80_cov["protein_id_short"].tolist()
    data_80_cov.to_csv(Path+fd[1]+"/data_80_cov.csv")
    if not os.path.exists(Path+fd[1]+"/80_cov_struct"):
        os.makedirs(Path+fd[1]+"/80_cov_struct")
    pdbl.download_pdb_files(list_80_cov, pdir=Path+fd[1]+"/80_cov_struct", obsolete=False)

    data_70_lig=fd[0].drop(fd[0][(fd[0]['coverage'] > 0.8)].index)
    data_70_lig=data_70_lig.drop(data_70_lig[(data_70_lig['Ligand_Similarity'] < 0.7)].index)
    list_70_lig=data_70_lig["protein_id_short"].tolist()
    data_70_lig.to_csv(Path+fd[1]+"/data_70_sim.csv")
    if not os.path.exists(Path+fd[1]+"/70_sim_struct"):
        os.makedirs(Path+fd[1]+"/70_sim_struct")
    pdbl.download_pdb_files(list_70_lig, pdir=Path+fd[1]+"/70_sim_struct", obsolete=False)
