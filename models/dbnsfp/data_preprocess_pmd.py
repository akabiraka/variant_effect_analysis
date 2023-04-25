import sys
sys.path.append("../variant_effect_analysis")
home_dir=""

import pandas as pd

from models.aa_common.data_loader import get_pmd_dataset


pmd_data_df = get_pmd_dataset(home_dir=home_dir)

pmd_data_df = pmd_data_df[~pd.isna(pmd_data_df["chrom_pos"])] # there are some rs-ids that does not have any data in dbSNP
print(pmd_data_df.shape)
# pmd_data_df = pmd_data_df[pd.isna(pmd_data_df["chrom_pos"])]
# print(pmd_data_df.shape)

pmd_data_df["chrom_pos"] = pmd_data_df["chrom_pos"].apply(lambda x: str(int(x))) # taking only chromosom number for dbNSFP inputs

print("Log: saving chromosomal variants PMD data ...")

pmd_data_df[["chrom", "chrom_pos", "ref_allele", "alt_allele"]].to_csv(f"models/dbnsfp/datasets_pmd/PMD_chromosomal.txt", index=False, sep=" ", header=False)


