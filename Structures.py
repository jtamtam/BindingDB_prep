import pandas as pd
from Bio.PDB import PDBList
from os import path
import os
from tqdm import tqdm


data_with_ligands_sorted=pd.read_csv("Fasta_files/data_with_ligands_sorted.csv",index_col=0)

#data = data.drop(data[(data['coverage'] != 1) & (data['Ligand_Similarity'] < 0.7)].index)

data_100_cov = data_with_ligands_sorted.drop(data_with_ligands_sorted[(data_with_ligands_sorted['coverage'] != 1)].index)
list_100_cov=data_100_cov["protein_id_short"].tolist()
data_100_cov.to_csv("data_100_cov.csv")

data_80_cov_70_sim = data_with_ligands_sorted.drop(data_with_ligands_sorted[(data_with_ligands_sorted['coverage'] < 0.8)].index)
data_80_cov_70_sim = data_80_cov_70_sim.drop(data_80_cov_70_sim[(data_80_cov_70_sim["Ligand_Similarity"] < 0.7)].index)
data_80_cov_70_sim = data_80_cov_70_sim.drop(data_80_cov_70_sim[(data_80_cov_70_sim["coverage"] == 1)].index)
list_80_cov_70_sim=data_80_cov_70_sim["protein_id_short"].tolist()
data_80_cov_70_sim.to_csv("data_80_cov_70_sim.csv")

data_80_cov = data_with_ligands_sorted.drop(data_with_ligands_sorted[(data_with_ligands_sorted['coverage'] < 0.8)].index)
data_80_cov = data_80_cov.drop(data_80_cov[(data_80_cov['Ligand_Similarity'] > 0.7)].index)
data_80_cov = data_80_cov.drop(data_80_cov[(data_80_cov['coverage'] == 1)].index)
list_80_cov=data_80_cov["protein_id_short"].tolist()
data_80_cov.to_csv("data_80_cov.csv")

data_70_lig=data_with_ligands_sorted.drop(data_with_ligands_sorted[(data_with_ligands_sorted['coverage'] > 0.8)].index)
data_70_lig=data_70_lig.drop(data_70_lig[(data_70_lig['Ligand_Similarity'] < 0.7)].index)
list_70_lig=data_70_lig["protein_id_short"].tolist()
data_70_lig.to_csv("data_70_lig.csv")

list_=data_with_ligands_sorted["protein_id_short"].tolist()

pdbl = PDBList()
if path.exists("Final_Output")==False: os.makedirs("Final_Structures")

pdbl.download_pdb_files(list_, pdir="Final_Output", obsolete=False)
