import sys
sys.path.append("../variant_effect_analysis")

import pandas as pd
import fileinput
from utils.ncbi_proteins import download_protein_list_mpi


def get_protein_mutation_df(inp_filepath):
    col_names = ["GRCh38Chromosome", "GRCh38Location", "Canonical SPDI"]
    df = pd.read_csv(inp_filepath, delim_whitespace=False, sep="\t")
    print(df.shape)


    prot_variations_df = df[["protein_accession.version", "Protein change"]]
    variations = []
    for i, row in enumerate(prot_variations_df.itertuples(index=False)):
        prot_acc_version, variants = row
        # print(prot_acc_version, variants)
        variants = [x.strip() for x in variants.split(",")]
        
        for v in variants:
            if v.endswith("fs"): 
                pos = int(v[1:-2])
                mut = v[-2:]
            elif v.endswith("del"): 
                pos = int(v[1:-3])
                mut = v[-3:]
            else: 
                pos = int(v[1:-1])
                mut = v[-1:]
                
                
            new_v = {"protein_accession.version": prot_acc_version,
                    "position": pos,
                    "wt": v[0],
                    "mut": mut}
            variations.append(new_v)
        # if i==10: break

    return pd.DataFrame(variations)



if __name__ == '__main__':
    inp_filepath = "data/clinvar/filtered/clinvar_HumanPathogenicMissenseVariants01012022To14022023.txt"
    # inp_filepath = "data/clinvar/filtered/clinvar_HumanLikelyPathogenicMissenseVariants01012022To14022023.txt"
    out_filepath = "models/aa_common/datasets_pathogenicity/" + inp_filepath.split("/")[-1].split(".")[0] + "_proteomicSNVs"
    
    prot_variations_df = get_protein_mutation_df(inp_filepath)
    prot_variations_df.to_csv(out_filepath+".txt", index=False, sep="\t", header=True)

    protein_acc_list = list(prot_variations_df["protein_accession.version"].unique())
    download_protein_list_mpi(protein_acc_list, len(protein_acc_list))
    print("#-unique NCBI protein sequences downloaded: ", len(protein_acc_list))
    
    
    print("Creating merged fasta document.")
    file_list = [f"data/proteins/fastas/{prot}.fasta" for prot in protein_acc_list]
    with open(out_filepath+".fasta", 'w') as file:
        input_lines = fileinput.input(file_list)
        file.writelines(input_lines)