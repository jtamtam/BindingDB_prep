from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from ssbio.protein.sequence.utils.blast import run_bidirectional_blast
from tqdm import tqdm
import os
import logging
logging.basicConfig(level = logging.INFO)

path="Fasta_files"
pdb_path=path+"/pdb/pdb.fasta"
if not os.path.exists("Fasta_files/BLAST_Output"):
    os.makedirs("Fasta_files/BLAST_Output")

for target in tqdm(os.listdir("Fasta_files/BindingDB")):
     filename = os.fsdecode(target)
     if filename.endswith(".fasta"):run_bidirectional_blast(path+"/BindingDB/"+filename,pdb_path, dbtype="prot",outdir="Fasta_files/BLAST_Output")
     else:continue