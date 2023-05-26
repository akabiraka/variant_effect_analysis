import sys
sys.path.append("../variant_effect_analysis")

import os
import pandas as pd
from Bio import SeqIO
# use only very common libraries here

def get_protein_sequences(home_dir="", max_seq_len=1022, return_type=None, data_type=None):
    """return_type: seq_record_list, protid_seq_tuple_list
       data_type: pathogenic, likely_pathogenic, popu_freq, pmd
    """
    print("\nLog: Loading combined fasta iterator ...")

    if data_type == "patho_and_likelypatho":
        filepath = home_dir+"models/aa_common/datasets_pathogenicity/patho_and_likelypatho.fasta"
    elif data_type == "pathogenic":
        filepath = home_dir+"models/aa_common/datasets_pathogenicity/clinvar_HumanPathogenicMissenseVariants01012022To14022023.fasta"
    elif data_type == "likely_pathogenic":
        filepath = home_dir+"models/aa_common/datasets_pathogenicity/clinvar_HumanLikelyPathogenicMissenseVariants01012022To14022023.fasta"
    elif data_type == "popu_freq":
        filepath = home_dir+"models/aa_common/datasets_population_freq/SNVs_with_popu_freq.fasta"
    elif data_type == "pmd":
        filepath = home_dir+"models/aa_common/datasets_pmd_analysis/pmd_sequences.fasta"
    else:
        raise NotImplementedError("'data_type' must be of pathogenic, likely_pathogenic, popu_freq, pmd")
        
        
    fasta_iterator = SeqIO.parse(filepath, format="fasta")
    
    if return_type == "seq_record_list":
        data = [seq_record for seq_record in fasta_iterator if len(str(seq_record.seq))<=max_seq_len]
    elif return_type == "protid_seq_tuple_list":
        data = [(seq_record.id, str(seq_record.seq)) for seq_record in fasta_iterator if len(str(seq_record.seq))<=max_seq_len]
    elif return_type == "protid_seq_dict":
        data = {seq_record.id: str(seq_record.seq) for seq_record in fasta_iterator if len(str(seq_record.seq))<=max_seq_len}
    else:
        raise NotImplementedError("'return_type' must be of seq_record_list, protid_seq_tuple_list")
    
    print(f"#-protein sequences (seq-len<={max_seq_len}): {len(data)}")
    return data
# get_protein_sequences()
# data = get_protein_sequences(home_dir="", max_seq_len=1022, return_type="protid_seq_tuple_list", data_type="pmd_analysis")
# print(data[0])


# ------------------------------the next 3 functions are for loading base 3 datasets-----------------------------
def get_pmd_dataset(home_dir=""):
    print("\nLog: Loading Protein Mutation Dataset (PMD) ...")
    pmd_df = pd.read_csv(home_dir+"models/aa_common/datasets_pmd_analysis/pmd.tsv", sep="\t") # PMD: protein mutation dataset
    pmd_df.drop_duplicates(keep="first", inplace=True, ignore_index=True)
    
    print(pmd_df.columns)
    print(pmd_df["functional_effect"].value_counts())
    print(pmd_df.shape)
    
    return pmd_df
# get_pmd_dataset()

def get_population_freq_SNVs(home_dir=""):
    print("\nLog: Loading data ...")
    variants_df = pd.read_csv(home_dir+"models/aa_common/datasets_population_freq/popu_freq.tsv", sep="\t")
    variants_df.drop_duplicates(keep="first", inplace=True, ignore_index=True)
    print(variants_df.columns)
    print(variants_df["class"].value_counts())
    print("total: ", variants_df.shape)

    return variants_df
# get_population_freq_SNVs()#force=True)


def get_patho_and_likelypatho_SNVs(home_dir=""):
    filepath = home_dir+f"models/aa_common/datasets_pathogenicity/patho_and_likelypatho.tsv"

    variants_df = pd.read_csv(filepath, sep="\t")
    variants_df.drop_duplicates(keep="first", inplace=True, ignore_index=True)
    print(variants_df.columns)
    print(variants_df["class"].value_counts())
    print("total: ", variants_df.shape)
    
    return variants_df
# get_patho_and_likelypatho_SNVs()


# the following 3 data-loader we are going to use for model running
def get_pmd_dbnsfp_dataset(home_dir=""):
    df = pd.read_csv(home_dir+f"models/aa_common/datasets_pmd_analysis/pmd_dbnsfp.tsv", sep="\t")

    fasta_iterator = SeqIO.parse(home_dir+f"models/aa_common/datasets_pmd_analysis/pmd_dbnsfp.fasta", format="fasta")
    protid_seq_dict = {seq_record.id: str(seq_record.seq) for seq_record in fasta_iterator}

    print(df.columns)
    print(df.shape)
    print(df["class"].value_counts())
    print("#-unique prots: ", len(protid_seq_dict))
    return df, protid_seq_dict
# get_pmd_dbnsfp_dataset()

def get_patho_likelypatho_neutral_dbnsfp_dataset(home_dir=""):
    pass

def get_popu_freq_dbnsfp_dataset(home_dir=""):
    pass