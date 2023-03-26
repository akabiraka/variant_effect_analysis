import sys
sys.path.append("../variant_effect_analysis")

import time
import numpy as np
import pandas as pd
from Bio import Entrez, SeqIO
from urllib.error import HTTPError
from mpi4py.futures import MPIPoolExecutor

from utils.ncbi_proteins import download_protein_list_mpi
home_dir=""#"../"

def get_prot_accessing_list(filepath, col_name):
    df = pd.read_csv(filepath, compression='gzip', delim_whitespace=False, sep="\t")#, header=None, names=col_names)
    protein_acc_list = list(df[col_name].unique())
    print("#-unique NCBI protein sequences to download: ", len(protein_acc_list))
    return protein_acc_list

if __name__ == '__main__':
    filepath = home_dir+"data/refseq/MANE.GRCh38.v1.0.summary.txt.gz" 
    col_name="RefSeq_prot"
    protein_acc_list = get_prot_accessing_list(filepath, col_name)
    
    # filepath = home_dir+"data/gene/gene2refseq_filtered.gz"
    # col_name="protein_accession.version"
    # protein_acc_list = get_prot_accessing_list(filepath, col_name)
    
    # protein_acc_list = protein_acc_list[:10]
    download_protein_list_mpi(protein_acc_list, len(protein_acc_list))

    # list_of_protein_acc_list = np.array_split(protein_acc_list, 1)
    
    # executor = MPIPoolExecutor()
    # # executor.map(download_prot_list, [protein_acc_list])
    # executor.map(download_prot_list, list_of_protein_acc_list)
    # # executor.map(download_protein_seq_and_save, protein_acc_list)
    
    # executor.shutdown()
     

# the following exact line of commands works
# salloc --partition=normal --mem=2G --ntasks=51
# module load gnu10 #gnu10/10.3.0-ya
# module load openmpi # openmpi/4.1.2-4a    
# source /projects/ashehu/akabir4/venvs/hopper_variant_effect_analysis_mine/bin/activate
# mpirun -np 21 python -m mpi4py.futures mpi/download_proteins.py