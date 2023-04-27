import sys
sys.path.append("../variant_effect_analysis")
home_dir = ""

import os
import numpy as np
import pandas as pd
from models.aa_common.data_loader import get_pmd_dataset

pmd_variants_df = get_pmd_dataset(home_dir=home_dir)
task = "pmd"

def filter_dbnsfp_preds(x):
    if x==".": return False
    elif ".;" in str(x): return False
    else: return True 
    
def compute_avg(x):
    x = str(x).split(";")
    return np.mean([float(i) for i in x if i!="."])    
    
def create_output_directories(model_name=None, task=None, home_dir=""):
    print("\nLog: Creating output directories ...") #-------------------------------------------
    model_out_dir = home_dir+f"models/dbnsfp/outputs/{model_name}/"
    model_task_out_dir = f"{model_out_dir}{task}/"
    os.makedirs(model_out_dir, exist_ok=True)
    os.makedirs(model_task_out_dir, exist_ok=True)
    return model_task_out_dir

def separate_dbnsfp_outputs_and_save(result_df, model_name, col_name):
    model_task_out_dir = create_output_directories(model_name, task, home_dir)
    # model_scores_df = result_df[result_df[col_name].apply(filter_dbnsfp_preds)] # removing those comparisons that does not produce any result
    model_scores_df = result_df.copy(deep=True)
    
    print(f"Log: Saving {model_name} prediction scores for {model_scores_df.shape[0]} SNVs...")
    model_scores = model_scores_df[col_name].apply(compute_avg) #lambda x: float(str(x).split(";")[0])) # can have multiple scores, ie '0.4573521;0.4573521;0.4573521;0.4573521'. taking 1st one
    result_df["pred"] = model_scores
    columns = ['mut_id', 'pmd_id', 'nr', 'crossref', 'uniprot_id', 'ensembl_id',
           'taxid', 'protein', 'mut_PMD', 'mut_real', 'wt', 'mut', 'prot_pos',
           'function_summarized', 'functional_effect', 'function', 'seq', 'snp_id',
           'mrna_acc', 'mrna_ver', 'mrna_pos', 'allele', 'protein_acc', 
           'protein_ver', 'verified', 'chrom', 'chrom_pos', 'variation',
           'variant_type', 'ref_allele', 'alt_allele', 'pmd_nr_id', 'pred']
    result_df = result_df[columns]
    result_df.to_csv(f"{model_task_out_dir}/preds_{model_name}.tsv", sep="\t", index=False, header=True)
    
    missing, total = result_df[pd.isna(result_df["pred"])].shape[0], result_df.shape[0]
    missing_values_percentage = (missing / total) * 100
    print(f"\tMissing values: ({missing}/{total})*100 = {missing_values_percentage}")
    print(f"Log: Saved df shape {model_name}-{col_name}: {result_df.shape}")
    


pred_df = pd.read_csv(home_dir+"models/dbnsfp/outputs/pmd_preds.txt", sep="\t")
pred_df = pred_df.loc[pred_df[["#chr", "pos(1-based)", "ref", "alt"]].drop_duplicates(keep="first").index] # for a single chromosomal position, a model can have multiple outputs from dbnsfp, so removing them
print(f"\n #-of SNVs found from dbNSFP: {pred_df.shape[0]}")


# pmd_variants_df["chrom_pos"] = pmd_variants_df[~pd.isna(pmd_variants_df["chrom_pos"])]["chrom_pos"].apply(lambda x: str(int(x))) # taking only chromosom number for dbNSFP inputs

print("Log: merging dbNSFP prediction scores with ground truth population freq data ...")
result_df = pmd_variants_df.merge(pred_df, how="left", left_on=["chrom", "chrom_pos", "ref_allele", "alt_allele"], right_on=["#chr", "pos(1-based)", "ref", "alt"])
result_df = result_df.drop_duplicates(keep="first")
print(result_df.columns)
print(result_df.shape)

separate_dbnsfp_outputs_and_save(result_df, model_name="sift", col_name="SIFT_score")
separate_dbnsfp_outputs_and_save(result_df, model_name="metarnn", col_name="MetaRNN_score")
separate_dbnsfp_outputs_and_save(result_df, model_name="mvp", col_name="MVP_score")
separate_dbnsfp_outputs_and_save(result_df, model_name="cadd", col_name="CADD_raw")
separate_dbnsfp_outputs_and_save(result_df, model_name="polyphen2_HVAR", col_name="Polyphen2_HVAR_score")
separate_dbnsfp_outputs_and_save(result_df, model_name="polyphen2_HDIV", col_name="Polyphen2_HDIV_score")
separate_dbnsfp_outputs_and_save(result_df, model_name="revel", col_name="REVEL_score")

