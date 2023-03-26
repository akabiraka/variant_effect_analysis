import sys
home_dir = ""
sys.path.append("../variant_effect_analysis")


import fileinput
import pandas as pd
from utils.ncbi_proteins import download_protein_list_mpi


def create_fasta(protein_acc_list):
    # the main input of this model is a fasta formatted sequences file. So the input file must only contain the relevant sequences.
    out_filepath = home_dir+"models/aa_common/datasets_population_freq/proteomic_SNVs"
    
    download_protein_list_mpi(protein_acc_list, len(protein_acc_list)) # download proteins
    print("#-unique NCBI protein sequences downloaded: ", len(protein_acc_list))

    print("Creating merged fasta document.")
    file_list = [home_dir+f"data/proteins/fastas/{prot}.fasta" for prot in protein_acc_list]
    with open(out_filepath+".fasta", 'w') as file:
        input_lines = fileinput.input(file_list)
        file.writelines(input_lines)


# TODO: check it again if it works
variants = pd.read_csv(home_dir+"models/aa_common/datasets_population_freq/proteomic_SNVs.txt", sep="\t")                
protein_acc_list = list(variants["prot_acc_version"].unique())
create_fasta(protein_acc_list)        