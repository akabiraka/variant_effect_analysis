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


def get_pmd_dataset(home_dir=""):
    print("\nLog: Loading Protein Mutation Dataset (PMD) ...")
    pmd_df = pd.read_csv(home_dir+"models/aa_common/datasets_pmd_analysis/pmd.tsv", sep="\t") # PMD: protein mutation dataset
    pmd_df.drop_duplicates(keep="first", inplace=True, ignore_index=True)
    
    # print("\nLog: excluding variants corresponding to proteins having seq-len>1022 ...")
    # protid_seq_tuple_list = get_protein_sequences(home_dir=home_dir, max_seq_len=1022, return_type="protid_seq_tuple_list", data_type="pmd")
    # new_protein_acc_list = list(zip(*protid_seq_tuple_list))[0]
    # pmd_df = pmd_df[pmd_df["pmd_nr_id"].isin(new_protein_acc_list)]
    
    print(pmd_df.columns)
    print(pmd_df["functional_effect"].value_counts())
    print(pmd_df.shape)
    
    return pmd_df
# get_pmd_dataset()

def get_population_freq_SNVs(home_dir=""):
    print("\nLog: Loading data ...")
    variants_df = pd.read_csv(home_dir+"models/aa_common/datasets_population_freq/SNVs_with_popu_freq.tsv", sep="\t")
    variants_df.drop_duplicates(keep="first", inplace=True, ignore_index=True)
    print(variants_df.columns)
    print(variants_df["class"].value_counts())
    print("total: ", variants_df.shape)

    return variants_df

    # data_filepath = home_dir+"models/aa_common/datasets_population_freq/SNVs_with_popu_freq_balanced.txt"
    # if os.path.exists(data_filepath) and not force:
    #     variants_df =  pd.read_csv(data_filepath, sep="\t")
    #     variants_df = variants_df.drop_duplicates(keep="first")
    #     print(variants_df.shape)
    #     print(variants_df.columns)
    #     print(variants_df["class"].value_counts())

    #     return variants_df
    

    # variants_df = pd.read_csv(home_dir+"models/aa_common/datasets_population_freq/SNVs_with_popu_freq.txt", sep="\t")
    # print(f"raw data: {variants_df.shape}")
    # print(variants_df.columns)
    # # print(variants_df.head())

    # print("\nLog: excluding variants corresponding to proteins having seq-len>1022 ...")
    # protid_seq_tuple_list = get_protein_sequences(home_dir="", max_seq_len=1022, return_type="protid_seq_tuple_list", data_type="popu_freq")
    # new_protein_acc_list = list(zip(*protid_seq_tuple_list))[0]
    # variants_df = variants_df[variants_df["prot_acc_version"].isin(new_protein_acc_list)]
    # print(variants_df.shape)

    # variants_df.loc[variants_df["mt_freq"]>=.01, "class"] = "Common"
    # variants_df.loc[(variants_df["mt_freq"]<.01) & (variants_df["mt_freq"]>=.001), "class"] = "Rare"
    # variants_df.loc[(variants_df["mt_freq"]<.001), "class"] = "Ultra-rare"
    # variants_df.loc[variants_df["mut_poulation"]==1, "class"] = "Singleton"
    # variants_df.loc[variants_df["mut_poulation"]<1, "class"] = "Zero-population"

    # filttered_variants_df = variants_df[variants_df["mut_poulation"]>1] # excluding singletons here

    # common = filttered_variants_df[filttered_variants_df["mt_freq"]>=.01]
    # rare = filttered_variants_df[(filttered_variants_df["mt_freq"]<.01) & (filttered_variants_df["mt_freq"]>=.001)]
    # ultra_rare = filttered_variants_df[filttered_variants_df["mt_freq"]<.001] # do not using this directly now
    # print(common.shape, rare.shape, ultra_rare.shape)


    # singletons = variants_df[variants_df["mut_poulation"]==1] # sampling singletons here
    # singletons_sampled = singletons.sample(n=common.shape[0]+rare.shape[0], random_state=1)
    # print(singletons_sampled.shape)

    # variants_df = pd.concat([common, rare, ultra_rare, singletons_sampled])
    # variants_df.reset_index(drop=True, inplace=True)
    # print(f"After combining common ({common.shape[0]}), rare ({rare.shape[0]}), ultra-rare ({ultra_rare.shape[0]}), sampled-singletons ({singletons.shape[0]}) and data: {variants_df.shape}")

    # variants_df = variants_df.drop_duplicates(keep="first")
    # variants_df.to_csv(data_filepath, sep="\t", header=True, index=False)

    # print(variants_df["class"].value_counts())
    # return variants_df

# get_population_freq_SNVs()#force=True)


def get_patho_and_likelypatho_SNVs(home_dir=""):
    filepath = home_dir+f"models/aa_common/datasets_pathogenicity/patho_and_likelypatho.tsv"

    variants_df = pd.read_csv(filepath, sep="\t")
    print(f"raw data: {variants_df.shape}")
    print(variants_df.columns)
    
    print("\nLog: excluding variants corresponding to proteins having seq-len>1022 ...")
    protid_seq_tuple_list = get_protein_sequences(home_dir=home_dir, max_seq_len=1022, return_type="protid_seq_tuple_list", data_type="patho_and_likelypatho")
    new_protein_acc_list = list(zip(*protid_seq_tuple_list))[0]
    variants_df = variants_df[variants_df["prot_acc_version"].isin(new_protein_acc_list)]
    
    variants_df = variants_df.drop_duplicates(keep="first")
    n_snp_ids = variants_df[~pd.isna(variants_df["snp_id"])].shape[0]
    print("#-of rs-ids mapped to pathogenicity dataset: ", n_snp_ids)
    
    print(variants_df["class"].value_counts())
    print(f"total: {variants_df.shape[0]}")
    return variants_df
# get_patho_and_likelypatho_SNVs()

# def map_NP_to_uniprot(df, col_name, home_dir=""):
#     np_to_uniprot_mapping_df = pd.read_csv(home_dir+"data/gene/np_to_uniprot_mapping.csv", sep="\t")
#     merged_df = pd.merge(left=df, right=np_to_uniprot_mapping_df, how="inner", left_on=col_name, right_on="NCBI_protein_accession")
#     merged_df.drop(columns="NCBI_protein_accession", inplace=True)
#     return merged_df

# def generate_neutral_SNVs(home_dir="", pathogenicity_type=None):
#     # pathogenicity_type: pathogenic, likely_pathogenic
#     variants_df = get_patho_and_likelypatho_SNVs(home_dir)
#     patho_variants_df = variants_df[variants_df["class"]==pathogenicity_type]

#     patho_unique_prot_acc_version_list = patho_variants_df["prot_acc_version"].unique()

#     print("\nLog: Loading population freq variants dataset ...")
#     popu_variants_df = get_population_freq_SNVs(home_dir)
#     print(popu_variants_df.columns)
    
#     commnon_popu_variants_df = popu_variants_df[popu_variants_df["class"]=="Common"]
#     commnon_popu_variants_df = commnon_popu_variants_df[commnon_popu_variants_df["prot_acc_version"].isin(patho_unique_prot_acc_version_list)] # variants must be of pathogenic (likely) proteins/genes
    
#     rare_popu_variants_df = popu_variants_df[popu_variants_df["class"]=="Rare"]
#     rare_popu_variants_df = rare_popu_variants_df[rare_popu_variants_df["prot_acc_version"].isin(patho_unique_prot_acc_version_list)] # variants must be of pathogenic (likely) proteins/genes
    
#     ultrarare_popu_variants_df = popu_variants_df[popu_variants_df["class"]=="Ultra-rare"]
#     ultrarare_popu_variants_df = ultrarare_popu_variants_df[ultrarare_popu_variants_df["prot_acc_version"].isin(patho_unique_prot_acc_version_list)] # variants must be of pathogenic (likely) proteins/genes
    
    
#     popu_variants_df = popu_variants_df[popu_variants_df["mut_poulation"]>1] # removing 0/1 population count variants
#     popu_variants_df = popu_variants_df[popu_variants_df["prot_acc_version"].isin(patho_unique_prot_acc_version_list)] # variants must be of pathogenic (likely) proteins/genes
#     print(f"{popu_variants_df.shape}")

#     out_dfs = []
#     for random_state in range(10):
#         data_filepath = home_dir+f"models/aa_common/datasets_pathogenicity/neutral_SNVs_for_{pathogenicity_type}_analysis/{random_state}.txt"
#         if os.path.exists(data_filepath):
#             sampled_neutral_variants_df = pd.read_csv(data_filepath, sep="\t")
#             print("Already exist", random_state, sampled_neutral_variants_df.shape)
#         else:
#             # sampled_neutral_variants_df = popu_variants_df.sample(n=patho_variants_df.shape[0], random_state=random_state)# sampling same number of pathogenic variants
#             # print(f"sampled population freq variants: {sampled_neutral_variants_df.shape}")
#             # sampled_neutral_variants_df["class"] = "neutral"

#             min(int(patho_variants_df.shape[0]/2), rare_popu_variants_df.shape[0])
#             sampled_commnon_popu_variants_df = commnon_popu_variants_df.sample(n=min(int(patho_variants_df.shape[0]/2), commnon_popu_variants_df.shape[0]), random_state=random_state)
#             sampled_rare_popu_variants_df = rare_popu_variants_df.sample(n=min(int(patho_variants_df.shape[0]/2), rare_popu_variants_df.shape[0]), random_state=random_state)
            
#             n_remaining = patho_variants_df.shape[0] - (sampled_commnon_popu_variants_df.shape[0] + sampled_rare_popu_variants_df.shape[0])
#             sampled_ultrarare_popu_variants_df = pd.DataFrame()
#             if n_remaining > 0:
#                 sampled_ultrarare_popu_variants_df = ultrarare_popu_variants_df.sample(n=n_remaining, random_state=random_state)
            
#             sampled_neutral_variants_df = pd.concat([sampled_commnon_popu_variants_df, sampled_rare_popu_variants_df, sampled_ultrarare_popu_variants_df])
#             sampled_neutral_variants_df = sampled_neutral_variants_df.rename(columns={"class": "popu_freq_class"})
#             sampled_neutral_variants_df["class"] = "neutral"

#             sampled_neutral_variants_df.to_csv(data_filepath, sep="\t", header=True, index=False)
#             print("Sampled", random_state, sampled_neutral_variants_df.shape)
            
#         out_dfs.append(sampled_neutral_variants_df)
#     return out_dfs

# generate_neutral_SNVs(home_dir="", pathogenicity_type="pathogenic")
# generate_neutral_SNVs(home_dir="", pathogenicity_type="likely_pathogenic")


# def get_pathogenicity_analysis_SNVs(home_dir="", pathogenicity_type=None):
#     """Deprecated. pathogenicity_type: pathogenic or likely_pathogenic, 
#     """
    
#     print("\nLog: Loading combined fasta iterator ...")
#     if pathogenicity_type == "pathogenic":
#         filepath = home_dir+"models/aa_common/datasets_pathogenicity/clinvar_HumanPathogenicMissenseVariants01012022To14022023.txt"
#     elif pathogenicity_type == "likely_pathogenic":
#         filepath = home_dir+"models/aa_common/datasets_pathogenicity/clinvar_HumanLikelyPathogenicMissenseVariants01012022To14022023.txt"
#     else:
#         raise NotImplementedError("'pathogenicity_type' must be of pathogenic, likely_pathogenic")
    
#     patho_variants_df = pd.read_csv(filepath, sep="\t")
#     print(f"{pathogenicity_type} raw data: {patho_variants_df.shape}")
#     print(patho_variants_df.columns)
    
    
#     print("\nLog: excluding variants corresponding to proteins having seq-len>1022 ...")
#     protid_seq_tuple_list = get_protein_sequences(home_dir=home_dir, max_seq_len=1022, return_type="protid_seq_tuple_list", data_type=pathogenicity_type)
#     new_protein_acc_list = list(zip(*protid_seq_tuple_list))[0]
#     patho_variants_df = patho_variants_df[patho_variants_df["prot_acc_version"].isin(new_protein_acc_list)]
#     print(patho_variants_df.shape)
    
#     patho_variants_df["id"] = patho_variants_df["clinvar_id"].apply(lambda x: "clinvar_id:"+str(x))
#     patho_variants_df = patho_variants_df[['id', 'chrom_acc_version', 'chrom_pos', 'ref_allele', 'alt_allele', 'prot_acc_version', 'prot_pos', 'wt', 'mut', "class"]]
    
#     patho_unique_prot_acc_version_list = patho_variants_df["prot_acc_version"].unique()
    
#     print("\nLog: Loading population freq variants dataset ...")
#     popu_variants_df = pd.read_csv(home_dir+"models/aa_common/datasets_population_freq/SNVs_with_popu_freq.txt", sep="\t")
#     print(popu_variants_df.columns)
    
#     popu_variants_df = popu_variants_df[popu_variants_df["mut_poulation"]>1] # removing 0/1 population count variants
#     popu_variants_df = popu_variants_df[popu_variants_df["prot_acc_version"].isin(patho_unique_prot_acc_version_list)] # variants must be of pathogenic (likely) proteins/genes
#     # print(f"{popu_variants_df.shape}")
    
    
#     out_dfs = []
#     for random_state in range(10):
#         data_filepath = home_dir+f"models/aa_common/datasets_pathogenicity/{pathogenicity_type}_and_neutral_SNVs/{random_state}.txt"
#         if os.path.exists(data_filepath):
#             variants_df = pd.read_csv(data_filepath, sep="\t")
#             print(random_state, variants_df.shape)
#         else:
#             sampled_popu_variants_df = popu_variants_df.sample(n=patho_variants_df.shape[0], random_state=random_state)# sampling same number of pathogenic variants
#             print(f"sampled population freq variants: {sampled_popu_variants_df.shape}")
#             sampled_popu_variants_df["id"] = sampled_popu_variants_df["snp_id"].apply(lambda x: "snp_id:"+str(x))
#             sampled_popu_variants_df["class"] = "neutral"
#             sampled_popu_variants_df = sampled_popu_variants_df[['id', 'chrom_acc_version', 'chrom_pos', 'ref_allele', 'alt_allele', 'prot_acc_version', 'prot_pos', 'wt', 'mut', "class"]]
            
#             variants_df = pd.concat([patho_variants_df, sampled_popu_variants_df], axis=0, ignore_index=True)
#             variants_df.reset_index(drop=True, inplace=True)
            
#             # print(variants_df["prot_acc_version"].value_counts())
#             variants_df.to_csv(data_filepath, sep="\t", header=True, index=False)
#             print(random_state, variants_df.shape)
            
#         out_dfs.append(variants_df)
#     return out_dfs

    
# get_pathogenicity_analysis_SNVs(pathogenicity_type="likely_pathogenic")
