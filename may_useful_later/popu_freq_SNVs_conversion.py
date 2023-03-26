import sys
sys.path.append("../variant_effect_analysis")    
home_dir = ""

import fileinput
import pandas as pd
import fileinput
from utils.ncbi_proteins import download_protein_list_mpi, download_protein_list, create_combined_fasta

from Bio import SeqIO
from Bio.PDB.Polypeptide import protein_letters_3to1 # 20 amino acids

def three_to_one(aa):
    if str.upper(aa) in protein_letters_3to1:
        return protein_letters_3to1[str.upper(aa)]
    return aa


def filter_unknown_variants(x):
    return len(x)==1



def get_variants_df(dbsnps_df):      
    # columns = ["snp_id", "chrom_acc_version", "chrom_pos", "ref_allele", "alt_allele", "prot_acc_version", "prot_pos", "wt", "mut", "wt_population", "mut_poulation", "wt_freq", "mt_freq"]
    variations = []
    for i, tuple in enumerate(dbsnps_df.itertuples()):
        # print(tuple.snp_id, tuple.variations, tuple.SAMN10492705)
        
        if len(tuple.REF)>1 or len(tuple.ALT)>1: # only considering single neucleodite variants
            continue
        
        prot_variations = tuple.variations.split(",") # ie: NP_064505.1:p.Arg898Lys,NP_064505.1:p.Arg898Met
        chrom_variations = tuple.ALT.split(",") # alt_alleles
        
        wt_population, mut_poulations = int(tuple.SAMN10492705.split(":")[0]), tuple.SAMN10492705.split(":")[1].split(",")
        total_population = wt_population+sum(list(map(int, mut_poulations)))
        
        
        try:
            for j, v in enumerate(prot_variations):
                if j < len(mut_poulations):
                    mut_poulation = int(mut_poulations[j])
                else: mut_poulation = 0

                if len(prot_variations) == len(chrom_variations):
                    # The chromosomal variants create same number of corresponding protein variants.
                    alt_allele = chrom_variations[j]
                else: alt_allele = chrom_variations[0]
                
                new_v = {"snp_id": tuple.snp_id,
                         
                         "chrom_acc_version": tuple.CHROM,
                         "chrom_pos": tuple.POS, # 1-indexed
                         "ref_allele": tuple.REF,
                         "alt_allele": alt_allele,

                         "prot_acc_version": v.split(":")[0], # protein_accession.version
                         "prot_pos": int(v.split(":")[1][5:-3]), # NCBI prot variants are 1-indexed
                         "wt": three_to_one(v.split(":")[1][2:5]),
                         "mut": three_to_one(v.split(":")[1][-3:]), 

                         "wt_population": wt_population,
                         "mut_poulation": mut_poulation, 
                         "wt_freq": wt_population/total_population, # freq should be computed here, since we are decomposing multiple SNVs into separate independent SNVs.
                         "mt_freq": mut_poulation/total_population}
                variations.append(new_v)
                # print(new_v)
        except:
            print(i, tuple)
            raise
        # if i==5000: break
    variations_df = pd.DataFrame(variations)
    return variations_df


if __name__ == "__main__":
    print("Log: loading raw variants data ...")
    inp_filepath = home_dir+"data/ALFA_population_freq/dbsnp_with_prots_and_population_freq_mapping.vcf.gz"
    dbsnps_df = pd.read_csv(inp_filepath, sep="\t", compression='gzip')
    # dbsnps_iterator = pd.read_csv(inp_filepath, sep="\t", compression='gzip', chunksize=2000)
    # dbsnps_df = dbsnps_iterator.__next__()
    print(dbsnps_df.shape)
    print(dbsnps_df.columns)
    # print(dbsnps_df.head())  


    print("\nLog: preprocessing raw variants data ...")
    variations_df = get_variants_df(dbsnps_df)
    print(variations_df.shape)


    print("\nLog: downloaing proteins ... ")
    protein_acc_list = list(variations_df["prot_acc_version"].unique())
    download_protein_list(protein_acc_list, start_i=0) # sequential downloading
    # download_protein_list_mpi(protein_acc_list, len(protein_acc_list))
    print("#-unique NCBI protein sequences downloaded: ", len(protein_acc_list))


    print("\nLog: excluding unknown amino acid variants ...")
    variations_df = variations_df[variations_df["wt"].apply(filter_unknown_variants)]
    variations_df = variations_df[variations_df["mut"].apply(filter_unknown_variants)]
    print(variations_df.shape)


    out_filepath = home_dir+"models/aa_common/datasets_population_freq/SNVs_with_popu_freq"
    print("\nLog: saving variants ...")
    variations_df.to_csv(out_filepath+".txt", index=False, sep="\t", header=True)
    

    print("\nLog: Creating merged fasta document ...")
    protein_acc_list = list(variations_df["prot_acc_version"].unique())
    create_combined_fasta(protein_acc_list, out_filepath+".fasta")