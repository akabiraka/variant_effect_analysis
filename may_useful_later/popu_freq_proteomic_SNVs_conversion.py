import sys
sys.path.append("../variant_effect_analysis")    
home_dir = ""

import fileinput
import pandas as pd
import fileinput
from utils.ncbi_proteins import download_protein_list_mpi, download_protein_list

from Bio import SeqIO
from Bio.PDB.Polypeptide import protein_letters_3to1 # 20 amino acids


def three_to_one(aa):
    if str.upper(aa) in protein_letters_3to1:
        return protein_letters_3to1[str.upper(aa)]
    return aa


def filter_unknown_variants(x):
    return len(x)==1


def create_fasta(protein_acc_list, out_filepath):
    file_list = [home_dir+f"data/proteins/fastas/{prot}.fasta" for prot in protein_acc_list]
    with open(out_filepath, 'w') as file:
        input_lines = fileinput.input(file_list)
        file.writelines(input_lines)
        
        
def get_protein_mutations_df(dbsnps_df):
    prot_variations = []
    for i, tuple in enumerate(dbsnps_df.itertuples()):
        # print(tuple.snp_id, tuple.variations, tuple.SAMN10492705)
        
        if len(tuple.REF)>1 or len(tuple.ALT)>1: # only considering single neucleodite variants
            continue
        
        variations = tuple.variations.split(",") # ie: NP_064505.1:p.Arg898Lys,NP_064505.1:p.Arg898Met
        
        wt_population, mut_poulations = int(tuple.SAMN10492705.split(":")[0]), tuple.SAMN10492705.split(":")[1].split(",")
        total_population = wt_population+sum(list(map(int, mut_poulations)))
        
        try:
            for j, v in enumerate(variations):
                if j < len(mut_poulations):
                    mut_poulation = int(mut_poulations[j])
                else: mut_poulation = 0
                
                new_v = {"prot_acc_version": v.split(":")[0], # protein_accession.version
                        "pos": int(v.split(":")[1][5:-3]),
                        "wt": three_to_one(v.split(":")[1][2:5]),
                        "mut": three_to_one(v.split(":")[1][-3:]), 
                        "wt_population": wt_population,
                        "mut_poulation": mut_poulation, 
                        "wt_freq": wt_population/total_population, # freq should be computed here, since we are decomposing multiple SNVs into separate independent SNVs.
                        "mt_freq": mut_poulation/total_population}
                prot_variations.append(new_v)
                # print(new_v)
        except:
            print(i, tuple)
            raise
        # if i==5000: break
    prot_variations_df = pd.DataFrame(prot_variations)
    return prot_variations_df




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
    prot_variations_df = get_protein_mutations_df(dbsnps_df)
    print(prot_variations_df.shape)


    print("\nLog: downloaing proteins ... ")
    protein_acc_list = list(prot_variations_df["prot_acc_version"].unique())
    download_protein_list(protein_acc_list, start_i=0) # sequential downloading
    # download_protein_list_mpi(protein_acc_list, len(protein_acc_list))
    print("#-unique NCBI protein sequences downloaded: ", len(protein_acc_list))
    
    
    print("\nLog: excluding proteins seq-len>1022 ...")
    new_protein_acc_list = []
    for prot_acc in protein_acc_list:
        seq_record = SeqIO.read(f"data/proteins/fastas/{prot_acc}.fasta", format="fasta") #
        if len(str(seq_record.seq)) <= 1022: 
            new_protein_acc_list.append(prot_acc)        
    print("#-of proteins after excluding (seq-len>1022): ", len(new_protein_acc_list))


    print("\nLog: excluding variants corresponding to proteins having seq-len>1022 ...")
    prot_variations_df = prot_variations_df[prot_variations_df["prot_acc_version"].isin(new_protein_acc_list)]
    print(prot_variations_df.shape)
    
    
    print("\nLog: excluding unknown amino acid variants ...")
    prot_variations_df = prot_variations_df[prot_variations_df["wt"].apply(filter_unknown_variants)]
    prot_variations_df = prot_variations_df[prot_variations_df["mut"].apply(filter_unknown_variants)]
    print(prot_variations_df.shape)
    
    
    out_filepath = home_dir+"models/aa_common/datasets_population_freq/proteomic_SNVs"
    print("\nLog: saving protein variants ...")
    prot_variations_df.to_csv(out_filepath+".txt", index=False, sep="\t", header=True)
    
    print("\nLog: Creating merged fasta document ...")
    protein_acc_list = list(prot_variations_df["prot_acc_version"].unique())
    create_fasta(protein_acc_list, out_filepath+".fasta")
        
        
    
# the following exact line of commands works
# salloc --partition=normal --mem=2G --ntasks=51
# module load gnu10 #gnu10/10.3.0-ya
# module load openmpi # openmpi/4.1.2-4a    
# source /projects/ashehu/akabir4/venvs/hopper_variant_effect_analysis_mine/bin/activate
# mpirun -np 41 python -m mpi4py.futures models/aa_common/proteomic_SNV_conversion_for_population_freq.py
# x must be <=41        