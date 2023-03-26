import sys
sys.path.append("../variant_effect_analysis")

import pandas as pd
from utils.pandas_extented_filters import filter_multi_base_variants, separate_ref_and_alt_allele

# inp_filepath = "data/clinvar/filtered/clinvar_HumanPathogenicMissenseVariants01012022To14022023.txt"
inp_filepath = "data/clinvar/filtered/clinvar_HumanLikelyPathogenicMissenseVariants01012022To14022023.txt"
col_names = ["GRCh38Chromosome", "GRCh38Location", "Canonical SPDI"]
df = pd.read_csv(inp_filepath, delim_whitespace=False, sep="\t")
print(df.shape)

chromosomal_variants_df = df[col_names]
temp = chromosomal_variants_df["Canonical SPDI"].apply(separate_ref_and_alt_allele) # filter 
chromosomal_variants_df["ref"], chromosomal_variants_df["alt"] = temp.apply(lambda x: x[0]), temp.apply(lambda x: x[1])

after_multibase_removal_df = chromosomal_variants_df[chromosomal_variants_df["ref"].apply(filter_multi_base_variants)] # filter: removing multi-nucleotide base variants 
after_multibase_removal_df = after_multibase_removal_df[after_multibase_removal_df["alt"].apply(filter_multi_base_variants)] # filter: removing multi-nucleotide base variants 
after_multibase_removal_df["GRCh38Chromosome"] = after_multibase_removal_df["GRCh38Chromosome"].apply(lambda x: x.split("|")[0])
after_multibase_removal_df

out_filepath = "models/aa_common/datasets_pathogenicity/" + inp_filepath.split("/")[-1].split(".")[0] + "_chromosomalSNVs.txt"
after_multibase_removal_df[["GRCh38Chromosome", "GRCh38Location", "ref", "alt"]].to_csv(out_filepath, index=False, sep=" ", header=False)