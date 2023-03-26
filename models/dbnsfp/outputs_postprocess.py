import sys
sys.path.append("../variant_effect_analysis")
home_dir = ""

import pandas as pd
from models.aa_common.data_loader import get_population_freq_SNVs

def filter_dbnsfp_preds(x):
    if x==".": return False
    elif ".;" in x: return False
    else: return True
    
variants_df = get_population_freq_SNVs(home_dir)    


pred_df = pd.read_csv(home_dir+"models/dbnsfp/outputs/chromosomal_SNVs_preds.txt", sep="\t")
print(f"#-of SNVs found from dbNSFP: {pred_df.shape[0]}")

variants_df["chrom"] = variants_df["chrom_acc_version"].apply(lambda x: int(x[x.index("_")+1:x.index(".")])) # taking only chromosom number for dbNSFP inputs

print("Log: merging dbNSFP prediction scores with ground truth population freq data ...")
result_df = pred_df.merge(variants_df, how="left", left_on=["#chr", "pos(1-based)", "ref", "alt"], right_on=["chrom", "chrom_pos", "ref_allele", "alt_allele"])
# result_df = variants_df.merge(pred_df, how="left", left_on=["chrom", "chrom_pos", "ref_allele", "alt_allele"], right_on=["#chr", "pos(1-based)", "ref", "alt"])
print(result_df.columns)


model_name = "metarnn"
metaRNN_scores_df = result_df[result_df["MetaRNN_score"].apply(filter_dbnsfp_preds)] # removing those comparisons that does not produce any result
print(f"Log: Saving MetaRNN prediction scores for {metaRNN_scores_df.shape[0]} SNVs...")
metaRNN_scores = metaRNN_scores_df["MetaRNN_score"].apply(lambda x: float(x.split(";")[0])) # can have multiple scores, ie '0.4573521;0.4573521;0.4573521;0.4573521'. taking 1st one
result_df["pred"] = metaRNN_scores
metaRNN_result_df = result_df[['snp_id', 'chrom_acc_version', 'chrom_pos', 'ref_allele', 'alt_allele', 'prot_acc_version', 'prot_pos', 'wt', 'mut', 'wt_population','mut_poulation', 'wt_freq', 'mt_freq', "pred"]]
metaRNN_result_df.to_csv(f"models/dbnsfp/outputs/popu_freq_preds_{model_name}.csv", index=False, sep="\t", header=True)


model_name = "mvp"
mvp_scores_df = result_df[result_df["MVP_score"].apply(filter_dbnsfp_preds)] # removing those comparisons that does not produce any result
print(f"Log: Saving MVP prediction scores for {mvp_scores_df.shape[0]} SNVs")
mvp_scores = mvp_scores_df["MVP_score"].apply(lambda x: float(x.split(";")[0]))
result_df["pred"] = mvp_scores
mvp_result_df = result_df[['snp_id', 'chrom_acc_version', 'chrom_pos', 'ref_allele', 'alt_allele', 'prot_acc_version', 'prot_pos', 'wt', 'mut', 'wt_population','mut_poulation', 'wt_freq', 'mt_freq', "pred"]]
mvp_result_df.to_csv(f"models/dbnsfp/outputs/popu_freq_preds_{model_name}.csv", index=False, sep="\t", header=True)