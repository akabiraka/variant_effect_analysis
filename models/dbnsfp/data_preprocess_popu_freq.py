import sys
sys.path.append("../variant_effect_analysis")
home_dir=""

import pandas as pd

from models.aa_common.data_loader import get_population_freq_SNVs


variants_df = get_population_freq_SNVs(home_dir=home_dir)

variants_df["chrom"] = variants_df["chrom_acc_version"].apply(lambda x: int(x[x.index("_")+1:x.index(".")])) # taking only chromosom number for dbNSFP inputs

print("Log: saving chromosomal variants w/o population freq ...")
variants_df[["chrom", "chrom_pos", "ref_allele", "alt_allele"]].to_csv(f"models/dbnsfp/datasets_population_freq/SNVs_with_popu_freq_balanced_1_chromosomal.txt", index=False, sep=" ", header=False)

