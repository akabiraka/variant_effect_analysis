import sys
sys.path.append("../variant_effect_analysis")

import pandas as pd
import fileinput
from utils.ncbi_proteins import download_protein_list_mpi, download_protein_list, create_combined_fasta

from Bio.PDB.Polypeptide import protein_letters_3to1 # 20 amino acids

def three_to_one(aa):
    if str.upper(aa) in protein_letters_3to1:
        return protein_letters_3to1[str.upper(aa)]
    return "unknown"

def get_variants_df(inp_filepath:str, pathogenicity_type=None) -> pd.DataFrame:   
    df = pd.read_csv(inp_filepath, sep="\t")
    # print(df.shape)
    # print(df.columns)
    # print(df.head())

    # pathogenicity_type: pathogenic, likely_pathogenic   
    # columns = ["Symbol", "GeneID", "GRCh38Chromosome", "GRCh38Location", "Canonical SPDI", "protein_accession.version", "Protein change"]

    variations = []
    for row_i in range(df.shape[0]):
        # print(df.loc[row_i])
        
        try:
            clinical_sig = df.loc[row_i, "Clinical significance (Last reviewed)"]
            if "Likely pathogenic" in clinical_sig:
                pathogenicity_type = "likely_pathogenic"
            else: pathogenicity_type = "pathogenic"
            
            # Case: bad->NM_032756.4(HPDL):c.[832G>A;91T>C], no protein variants, good->NM_000478.6(ALPL):c.1276G>A (p.Gly426Ser)
            x = df.loc[row_i, "Name"].split() # 
            if len(x)!=2: continue 
            
            prot_variant = x[1]
            try:
                prot_pos = int(prot_variant[6:-4]) # case: NM_003001.5(SDHC):c.452_455delinsATGA (p.Ser151_Gly152delinsTyrGlu)
            except:
                continue
            
            wt, mut = prot_variant[3:6], prot_variant[-4:-1]
            wt, mut = three_to_one(wt), three_to_one(mut)
            if wt=="unknown" or mut=="unknown": continue # not considering unknown protein variants

            x = df.loc[row_i, "Canonical SPDI"] 
            x = x.split("|")[0] # Case: NC_000023.11:154532268:C:A|NC_000023.11:154532045:A:C. Take the 1st.
            chrom_acc_version, _, ref_allele, alt_allele = x.split(":") # ie: NC_000023.11:154532268:C:A # chr-pos in Canonical SPDI is 0-indexed
            if len(ref_allele)>1 or len(alt_allele)>1: continue # only considering single neucleodite variants
            
            chrom_pos = int(df.loc[row_i, "GRCh38Location"]) #1-indexed
            
            
            
            new_v = {"clinvar_id": df.loc[row_i, "VariationID"], 
                    "gene_symbol": df.loc[row_i, "Symbol"],
                    "gene_id": df.loc[row_i, "GeneID"], 
                    
                    "chrom_acc_version": chrom_acc_version,
                    "chrom_pos": chrom_pos, # 1-indexed
                    "ref_allele": ref_allele,
                    "alt_allele": alt_allele,
                    
                    "prot_acc_version": df.loc[row_i, "protein_accession.version"],
                    "prot_pos": prot_pos, # NCBI prot variants are 1-indexed
                    "wt": wt, # 1-letter amino acid
                    "mut": mut,
                    "class": pathogenicity_type
                }
        
            variations.append(new_v) 
        except:
            print("Error occured: ")
            print(row_i, df.loc[row_i])
            raise
        
        # if row_i==10: break
    variations_df = pd.DataFrame(variations)
    return variations_df                      
            
        
        
if __name__ == '__main__':
    inp_filepath = "data/clinvar/filtered/clinvar_HumanPathogenicMissenseVariants01012022To14022023.txt"
    patho_variations_df = get_variants_df(inp_filepath, pathogenicity_type="pathogenic")
    print(patho_variations_df.shape)

    inp_filepath = "data/clinvar/filtered/clinvar_HumanLikelyPathogenicMissenseVariants01012022To14022023.txt"
    likelypatho_variations_df = get_variants_df(inp_filepath, pathogenicity_type="likely_pathogenic")
    print(likelypatho_variations_df.shape)
    
    variations_df = pd.concat([patho_variations_df, likelypatho_variations_df], ignore_index=True)
    variations_df = variations_df.drop_duplicates(keep="first")
    n_patho = variations_df[variations_df["class"] == "pathogenic"].shape[0]
    n_likelypatho = variations_df[variations_df["class"] == "likely_pathogenic"].shape[0]
    print(f"Pathogenic: {n_patho}, Likely-pathogenic: {n_likelypatho}")
    
    
    # print("\nLog: downloaing proteins ... ")
    # protein_acc_list = list(variations_df["prot_acc_version"].unique())
    # download_protein_list(protein_acc_list, start_i=0) # sequential downloading
    # # download_protein_list_mpi(protein_acc_list, len(protein_acc_list))
    # print("#-unique NCBI protein sequences downloaded: ", len(protein_acc_list))
    

    filename = "patho_and_likelypatho"
    out_filepath = f"models/aa_common/datasets_pathogenicity/{filename}"
    print("\nLog: saving variants ...")
    variations_df.to_csv(out_filepath+".txt", index=False, sep="\t", header=True)
    

    # print("\nLog: Creating merged fasta document ...")
    # protein_acc_list = list(variations_df["prot_acc_version"].unique())
    # create_combined_fasta(protein_acc_list, out_filepath+".fasta")
    print(variations_df.shape)












# ----------------------------------------------------------------
# prot_variants = df.loc[row_i, "Protein change"]
# prot_variants = [x.strip() for x in prot_variants.split(",")]
# try:   
#     for v in prot_variants:
#         if v.endswith("fs") or v.endswith("del"): continue # not considering framshift (fs) or deletion (del) variants
        
#         x = df.loc[row_i, "Canonical SPDI"] 
#         x = x.split("|")[0] # Case: NC_000023.11:154532268:C:A|NC_000023.11:154532045:A:C. Take the 1st.
#         chrom_acc_version, _, ref_allele, alt_allele = x.split(":") # ie: NC_000023.11:154532268:C:A # chr-pos in Canonical SPDI is 0-indexed
#         chrom_pos = int(df.loc[row_i, "GRCh38Location"]) #1-indexed
        
#         if len(ref_allele)>1 or len(alt_allele)>1: continue # only considering single neucleodite variants

#         new_v = {"gene_symbol": df.loc[row_i, "Symbol"],
#                 "gene_id": df.loc[row_i, "GeneID"], 
                
#                 "chrom_acc_version": chrom_acc_version,
#                 "chrom_pos": chrom_pos, # 1-indexed
#                 "ref_allele": ref_allele,
#                 "alt_allele": alt_allele,
                
#                 "prot_acc_version": df.loc[row_i, "protein_accession.version"],
#                 "prot_pos": int(v[1:-1]), # NCBI prot variants are 1-indexed
#                 "wt": v[0], # 1-letter amino acid
#                 "mut": v[-1:]}
        
#         variations.append(new_v)