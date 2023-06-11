import sys
sys.path.append("../variant_effect_analysis")

import pandas as pd
from Bio import SeqIO
# use only very common libraries here


def get_protein_sequences(task, return_type=None, home_dir=""):
    """
    task: pmd, patho, popu_freq
    seq_return_type: protid_seq_dict, seq_record_list
    """

    if task=="pmd":
        filepath = home_dir+f"models/aa_common/datasets_pmd_analysis/pmd_dbnsfp.fasta"
    if task=="patho":
        filepath = home_dir+f"models/aa_common/datasets_pathogenicity/patho_likelypatho_neutral_dbnsfp.fasta"
    elif task=="popu_freq":
        filepath = home_dir+f"models/aa_common/datasets_population_freq/popu_freq_with_dbnsfp_sampled.fasta"

    fasta_iterator = SeqIO.parse(filepath, format="fasta")

    if return_type == "seq_record_list":
        data = [seq_record for seq_record in fasta_iterator]
    else: # protid_seq_dict
        data = {seq_record.id: str(seq_record.seq) for seq_record in fasta_iterator}

    return data
# get_protein_sequences()



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


# ----------------------------------------the following 3 data-loader we are going to use for model running------------------------
def get_pmd_dbnsfp_dataset(home_dir="", seq_return_type=None):
    df = pd.read_csv(home_dir+f"models/aa_common/datasets_pmd_analysis/pmd_dbnsfp.tsv", sep="\t")
    seq_data = get_protein_sequences(task="pmd", return_type=seq_return_type, home_dir=home_dir)

    print(df.columns)
    print(df.shape)
    print(df["class"].value_counts())
    print("#-unique prots: ", len(seq_data))
    return df, seq_data
# get_pmd_dbnsfp_dataset(seq_return_type="seq_record_list")

def get_patho_likelypatho_neutral_dbnsfp_dataset(home_dir="", seq_return_type=None):
    df = pd.read_csv(home_dir+f"models/aa_common/datasets_pathogenicity/patho_likelypatho_neutral_dbnsfp.tsv", sep="\t")
    seq_data = get_protein_sequences(task="patho", return_type=seq_return_type, home_dir=home_dir)

    print(df.columns)
    print(df.shape)
    print(df["class"].value_counts())
    print("#-unique prots: ", len(seq_data))
    return df, seq_data
# get_patho_likelypatho_neutral_dbnsfp_dataset()

def get_popu_freq_dbnsfp_dataset(home_dir="", seq_return_type=None):
    df = pd.read_csv(home_dir+f"models/aa_common/datasets_population_freq/popu_freq_with_dbnsfp_sampled.tsv", sep="\t")
    seq_data = get_protein_sequences(task="popu_freq", return_type=seq_return_type, home_dir=home_dir)

    print(df.columns)
    print(df.shape)
    print(df["class"].value_counts())
    print("#-unique prots: ", len(seq_data))
    return df, seq_data
# get_popu_freq_dbnsfp_dataset()