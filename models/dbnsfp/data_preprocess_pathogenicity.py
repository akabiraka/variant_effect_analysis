import sys
sys.path.append("../variant_effect_analysis")
home_dir=""

import pandas as pd

from models.aa_common.data_loader import get_pathogenicity_analysis_SNVs

pathogenicity_type = "likely_pathogenic" # pathogenic, likely_pathogenic
list_of_variants_df = get_pathogenicity_analysis_SNVs(home_dir=home_dir, pathogenicity_type=pathogenicity_type)

for i, variants_df in enumerate(list_of_variants_df):
    variants_df["chrom"] = variants_df["chrom_acc_version"].apply(lambda x: int(x[x.index("_")+1:x.index(".")])) # taking only chromosom number for dbNSFP inputs

    print("Log: saving chromosomal variants w/o population freq ...")
    variants_df[["chrom", "chrom_pos", "ref_allele", "alt_allele"]].to_csv(f"models/dbnsfp/datasets_pathogenicity/{pathogenicity_type}_and_neutral_SNVs/{i}.txt", index=False, sep=" ", header=False)

    # break