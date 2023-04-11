import sys
sys.path.append("../variant_effect_analysis")
home_dir=""

import pandas as pd

from models.aa_common.data_loader import get_patho_and_likelypatho_SNVs

task = "patho_and_likelypatho"
patho_and_likelypatho_variants_df = get_patho_and_likelypatho_SNVs(home_dir)

patho_and_likelypatho_variants_df["chrom"] = patho_and_likelypatho_variants_df["chrom_acc_version"].apply(lambda x: int(x[x.index("_")+1:x.index(".")])) # taking only chromosom number for dbNSFP inputs

print("Log: saving chromosomal variants w/o population freq ...")
patho_and_likelypatho_variants_df[["chrom", "chrom_pos", "ref_allele", "alt_allele"]].to_csv(f"models/dbnsfp/datasets_pathogenicity/{task}.txt", index=False, sep=" ", header=False)