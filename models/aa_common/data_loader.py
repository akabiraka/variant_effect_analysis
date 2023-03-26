import sys
sys.path.append("../variant_effect_analysis")

import os
import pandas as pd
from Bio import SeqIO
# use only very common libraries here

def get_protein_sequences(home_dir="", max_seq_len=1022, return_type=None, data_type=None):
    """return_type: seq_record_list, protid_seq_tuple_list
       data_type: pathogenic, likely_pathogenic, popu_freq, pmd_analysis
    """
    print("\nLog: Loading combined fasta iterator ...")

    if data_type == "pathogenic":
        filepath = home_dir+"models/aa_common/datasets_pathogenicity/clinvar_HumanPathogenicMissenseVariants01012022To14022023.fasta"
    elif data_type == "likely_pathogenic":
        filepath = home_dir+"models/aa_common/datasets_pathogenicity/clinvar_HumanLikelyPathogenicMissenseVariants01012022To14022023.fasta"
    elif data_type == "popu_freq":
        filepath = home_dir+"models/aa_common/datasets_population_freq/SNVs_with_popu_freq.fasta"
    elif data_type == "pmd_analysis":
        filepath = home_dir+"models/aa_common/datasets_pmd_analysis/pmd_sequences.fasta"
    else:
        raise NotImplementedError("'data_type' must be of pathogenic, likely_pathogenic, popu_freq")
        
        
    fasta_iterator = SeqIO.parse(filepath, format="fasta")
    
    if return_type == "seq_record_list":
        data = [seq_record for seq_record in fasta_iterator if len(str(seq_record.seq))<=max_seq_len]
    elif return_type == "protid_seq_tuple_list":
        data = [(seq_record.id, str(seq_record.seq)) for seq_record in fasta_iterator if len(str(seq_record.seq))<=max_seq_len]
    else:
        raise NotImplementedError("'return_type' must be of seq_record_list, protid_seq_tuple_list")
    
    print(f"#-protein sequences (seq-len<={max_seq_len}): {len(data)}")
    return data
# get_protein_sequences()


def get_pmd_analysis_dataset(home_dir=""):
    print("\nLog: Loading Protein Mutation Dataset (PMD) ...")
    pmd_df = pd.read_csv(home_dir+"models/aa_common/datasets_pmd_analysis/pmd_data.csv", sep=",") # PMD: protein mutation dataset
    print(pmd_df.shape)
    print(pmd_df.columns)
    
    print("\nLog: excluding variants corresponding to proteins having seq-len>1022 ...")
    protid_seq_tuple_list = get_protein_sequences(home_dir=home_dir, max_seq_len=1022, return_type="protid_seq_tuple_list", data_type="pmd_analysis")
    new_protein_acc_list = list(zip(*protid_seq_tuple_list))[0]
    pmd_df = pmd_df[pmd_df["protein_id"].isin(new_protein_acc_list)]
    print(pmd_df.shape)
    
    return pmd_df
# get_pmd_analysis_dataset()


def get_pathogenicity_analysis_SNVs(home_dir="", pathogenicity_type=None):
    """pathogenicity_type: pathogenic or likely_pathogenic, 
    """
    
    print("\nLog: Loading combined fasta iterator ...")
    if pathogenicity_type == "pathogenic":
        filepath = home_dir+"models/aa_common/datasets_pathogenicity/clinvar_HumanPathogenicMissenseVariants01012022To14022023.txt"
    elif pathogenicity_type == "likely_pathogenic":
        filepath = home_dir+"models/aa_common/datasets_pathogenicity/clinvar_HumanLikelyPathogenicMissenseVariants01012022To14022023.txt"
    else:
        raise NotImplementedError("'pathogenicity_type' must be of pathogenic, likely_pathogenic")
    
    patho_variants_df = pd.read_csv(filepath, sep="\t")
    print(f"{pathogenicity_type} raw data: {patho_variants_df.shape}")
    print(patho_variants_df.columns)
    
    
    print("\nLog: excluding variants corresponding to proteins having seq-len>1022 ...")
    protid_seq_tuple_list = get_protein_sequences(home_dir=home_dir, max_seq_len=1022, return_type="protid_seq_tuple_list", data_type=pathogenicity_type)
    new_protein_acc_list = list(zip(*protid_seq_tuple_list))[0]
    patho_variants_df = patho_variants_df[patho_variants_df["prot_acc_version"].isin(new_protein_acc_list)]
    print(patho_variants_df.shape)
    
    patho_variants_df["id"] = patho_variants_df["clinvar_id"].apply(lambda x: "clinvar_id:"+str(x))
    patho_variants_df = patho_variants_df[['id', 'chrom_acc_version', 'chrom_pos', 'ref_allele', 'alt_allele', 'prot_acc_version', 'prot_pos', 'wt', 'mut', "class"]]
    
    patho_unique_prot_acc_version_list = patho_variants_df["prot_acc_version"].unique()
    
    print("\nLog: Loading population freq variants dataset ...")
    popu_variants_df = pd.read_csv(home_dir+"models/aa_common/datasets_population_freq/SNVs_with_popu_freq.txt", sep="\t")
    print(popu_variants_df.columns)
    
    popu_variants_df = popu_variants_df[popu_variants_df["mut_poulation"]>1] # removing 0/1 population count variants
    popu_variants_df = popu_variants_df[popu_variants_df["prot_acc_version"].isin(patho_unique_prot_acc_version_list)] # variants must be of pathogenic (likely) proteins/genes
    # print(f"{popu_variants_df.shape}")
    
    
    out_dfs = []
    for random_state in range(10):
        data_filepath = home_dir+f"models/aa_common/datasets_pathogenicity/{pathogenicity_type}_and_neutral_SNVs/{random_state}.txt"
        if os.path.exists(data_filepath):
            variants_df = pd.read_csv(data_filepath, sep="\t")
            print(random_state, variants_df.shape)
        else:
            sampled_popu_variants_df = popu_variants_df.sample(n=patho_variants_df.shape[0], random_state=random_state)# sampling same number of pathogenic variants
            print(f"sampled population freq variants: {sampled_popu_variants_df.shape}")
            sampled_popu_variants_df["id"] = sampled_popu_variants_df["snp_id"].apply(lambda x: "snp_id:"+str(x))
            sampled_popu_variants_df["class"] = "neutral"
            sampled_popu_variants_df = sampled_popu_variants_df[['id', 'chrom_acc_version', 'chrom_pos', 'ref_allele', 'alt_allele', 'prot_acc_version', 'prot_pos', 'wt', 'mut', "class"]]
            
            variants_df = pd.concat([patho_variants_df, sampled_popu_variants_df], axis=0, ignore_index=True)
            variants_df.reset_index(drop=True, inplace=True)
            
            # print(variants_df["prot_acc_version"].value_counts())
            variants_df.to_csv(data_filepath, sep="\t", header=True, index=False)
            print(random_state, variants_df.shape)
            
        out_dfs.append(variants_df)
    return out_dfs

    
# get_pathogenicity_analysis_SNVs(pathogenicity_type="likely_pathogenic")





def get_population_freq_SNVs(home_dir="", force=False):
    print("\nLog: Loading data ...")
    data_filepath = home_dir+"models/aa_common/datasets_population_freq/SNVs_with_popu_freq_balanced.txt"
    if os.path.exists(data_filepath) and not force:
        variants_df =  pd.read_csv(data_filepath, sep="\t")
        variants_df = variants_df.drop_duplicates(keep="first")
        print(variants_df.shape)
        print(variants_df.columns)

        common = variants_df[variants_df["mt_freq"]>=.01]
        rare = variants_df[(variants_df["mt_freq"]<.01) & (variants_df["mt_freq"]>=.001)]
        singletons = variants_df[variants_df["mut_poulation"]==1]
        print(f"After combining common ({common.shape[0]}), rare ({rare.shape[0]}) and sampled-singletons ({singletons.shape[0]}), data: {variants_df.shape}")
        print(variants_df["prot_acc_version"].value_counts())

        return variants_df
    

    variants_df = pd.read_csv(home_dir+"models/aa_common/datasets_population_freq/SNVs_with_popu_freq.txt", sep="\t")
    print(f"raw data: {variants_df.shape}")
    print(variants_df.columns)
    # print(variants_df.head())

    print("\nLog: excluding variants corresponding to proteins having seq-len>1022 ...")
    protid_seq_tuple_list = get_protein_sequences(home_dir="", max_seq_len=1022, return_type="protid_seq_tuple_list", data_type="popu_freq")
    new_protein_acc_list = list(zip(*protid_seq_tuple_list))[0]
    variants_df = variants_df[variants_df["prot_acc_version"].isin(new_protein_acc_list)]
    print(variants_df.shape)


    filttered_variants_df = variants_df[variants_df["mut_poulation"]>=1]
    common = filttered_variants_df[filttered_variants_df["mt_freq"]>=.01]
    rare = filttered_variants_df[(filttered_variants_df["mt_freq"]<.01) & (filttered_variants_df["mt_freq"]>=.001)]
    # ultra_rare = filttered_variants_df[filttered_variants_df["mt_freq"]<.001] # do not using this directly now
    # print(common.shape, rare.shape, ultra_rare.shape)


    singletons = variants_df[variants_df["mut_poulation"]==1]
    # print(singletons.shape)

    singletons_sampled = singletons.sample(n=common.shape[0]+rare.shape[0]) # random_state=1
    # print(singletons_sampled.shape)

    variants_df = pd.concat([common, rare, singletons_sampled])
    variants_df.reset_index(drop=True, inplace=True)
    print(f"After combining common ({common.shape[0]}), rare ({rare.shape[0]}) and sampled-singletons ({singletons_sampled.shape[0]}), data: {variants_df.shape}")

    # print(variants_df[variants_df["mt_freq"].isna()])
    print(variants_df["prot_acc_version"].value_counts())
    variants_df = variants_df.drop_duplicates(keep="first")
    variants_df.to_csv(data_filepath, sep="\t", header=True, index=False)

    return variants_df

# get_population_freq_SNVs()#force=True)
