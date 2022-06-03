import array
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import wget
import os
from tqdm import tqdm
from collections import defaultdict

if not os.path.exists("Final_Data"):
    os.makedirs("Final_Data")

if not os.path.exists("cc-to-pdb.tdd"):
    wget.download("http://ligand-expo.rcsb.org/dictionaries/cc-to-pdb.tdd", "cc-to-pdb.tdd")


read_file = pd.read_csv (r'cc-to-pdb.tdd',sep='\t', names=['Keys', 'Lists'])
read_file['Lists']=read_file['Lists'].str.rstrip(" ")
read_file['Lists']=read_file['Lists'].str.split(" ")
lig_to_prot = dict(zip(read_file.Keys, read_file.Lists))
prot_to_lig = {}
for k in lig_to_prot:
    for v in lig_to_prot[k]:
        prot_to_lig.setdefault(v, []).append(k)

def tanimoto_calc(smi1, smi2):
    mol1 = Chem.MolFromSmiles(smi1)
    mol2 = Chem.MolFromSmiles(smi2)
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 3, nBits=2048)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 3, nBits=2048)
    s = round(DataStructs.TanimotoSimilarity(fp1,fp2),3)
    return s

if not os.path.exists("Components-smiles-stereo-oe.smi"):
    wget.download("http://ligand-expo.rcsb.org/dictionaries/Components-smiles-stereo-oe.smi", "Components-smiles-stereo-oe.smi")

ligand_data=pd.read_csv("Components-smiles-stereo-oe.smi",sep="\t",header=None, names=["Chemical","Ligand_ID", "Ligand"],dtype=str)

ligand_data["Chemical"]=ligand_data["Chemical"].astype(str)
ligand_data["Ligand"]=ligand_data["Ligand"].astype(str)
ligand_data.set_index("Ligand_ID",drop=False,inplace=True)
ligand_data_dict=ligand_data.to_dict("index")


train_ic50=pd.read_csv("Before_Ligands/train_ic50.csv",index_col=0)
test_ic50=pd.read_csv("Before_Ligands/test_ic50.csv",index_col=0)
valid_ic50=pd.read_csv("Before_Ligands/valid_ic50.csv",index_col=0)

train_ki=pd.read_csv("Before_Ligands/train_ki.csv",index_col=0)
test_ki=pd.read_csv("Before_Ligands/test_ki.csv",index_col=0)
valid_ki=pd.read_csv("Before_Ligands/valid_ki.csv",index_col=0)

train_kd=pd.read_csv("Before_Ligands/train_kd.csv",index_col=0)
test_kd=pd.read_csv("Before_Ligands/test_kd.csv",index_col=0)
valid_kd=pd.read_csv("Before_Ligands/valid_kd.csv",index_col=0)

array_before_ligands=[
    train_ic50,
    test_ic50,
    valid_ic50,
    train_ki,
    test_ki,
    valid_ki,
    train_kd,
    test_kd,
    valid_kd,
    ]

for abl in array_before_ligands:
    abl["Ligand_Similarity"]=""
    abl["Ligand_ID"]=""
    abl["Ligand_Formula"]=""
    abl.reset_index(drop=True, inplace=True) 
    
array_before_ligand_dict=[]
for abl in array_before_ligands:
    array_before_ligand_dict.append(abl.to_dict("index"))

#%%
array_after_ligand=[]
for abl in tqdm(array_before_ligand_dict):
    print("longeur=",len(abl))
    for k, v in tqdm(abl.items()):
        if prot_to_lig.get(abl[k]["protein_id_short"]) is None:

            abl[k]["Ligand_Similarity"]=float(0)
            abl[k]["Ligand_ID"]="/"
            abl[k]["Ligand_Formula"]="/"
            continue 
        else:
            ligand_to_compare=prot_to_lig.get(abl[k]["protein_id_short"])
            max_tani=-100000000
            saved_row=None
            for ele in ligand_to_compare:
                
                if ligand_data_dict.get(ele) is None:
                    abl[k]["Ligand_Similarity"]=float(0)
                    abl[k]["Ligand_ID"]="/"
                    abl[k]["Ligand_Formula"]="/"
                    continue
                row=ligand_data_dict[ele]
               
                lig_from_pdb=row["Chemical"]
                lig_from_data=abl[k]["Drug"]
                try:
                    maybe_max=tanimoto_calc(lig_from_pdb,lig_from_data)
                except:
                    pass
                if maybe_max>max_tani: 
                    max_tani=maybe_max
                    saved_row=row
            abl[k]["Ligand_Similarity"]=float(max_tani)
            abl[k]["Ligand_Formula"]=saved_row["Chemical"]
            abl[k]["Ligand_ID"]=saved_row["Ligand_ID"]

    abl_df=pd.DataFrame.from_dict(abl,orient="index")
    abl_df = abl_df.drop(abl_df[(abl_df['coverage'] < 0.8) & (abl_df['Ligand_Similarity'] < 0.7)].index)
    abl_df.reset_index(drop=True, inplace=True)
    array_after_ligand.append(abl_df)



array_after_ligand[0].to_csv("Final_Data/train_ic50_final.csv")
array_after_ligand[1].to_csv("Final_Data/test_ic50_final.csv")
array_after_ligand[2].to_csv("Final_Data/valid_ic50_final.csv")

array_after_ligand[3].to_csv("Final_Data/train_ki_final.csv")
array_after_ligand[4].to_csv("Final_Data/test_ki_final.csv")
array_after_ligand[5].to_csv("Final_Data/valid_ki_final.csv")

array_after_ligand[6].to_csv("Final_Data/train_kd_final.csv")
array_after_ligand[7].to_csv("Final_Data/test_kd_final.csv")
array_after_ligand[8].to_csv("Final_Data/valid_kd_final.csv")
