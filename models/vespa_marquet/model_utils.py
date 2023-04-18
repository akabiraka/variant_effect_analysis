import sys
sys.path.append("../variant_effect_analysis")

import os
import numpy as np

import utils.pickle_utils as pickle_utils
import h5py




def create_output_directories(model_name=None, task=None, home_dir=""):
    """model_name: protbert, unirep
       task: pathogenic, likely_pathogenic
    """
    print("\nLog: Creating output directories ...") #-------------------------------------------
    model_out_dir = home_dir+f"models/vespa_marquet/outputs/{model_name}/"
    model_logits_out_dir = f"{model_out_dir}lm_outputs/"
    model_task_out_dir = f"{model_out_dir}{task}/"
    os.makedirs(model_out_dir, exist_ok=True)
    os.makedirs(model_logits_out_dir, exist_ok=True)
    os.makedirs(model_task_out_dir, exist_ok=True)
    return model_task_out_dir, model_logits_out_dir

filepath = "models/vespa_marquet/cache/vespal_human_proteome.h5"
computed_data_file_handle = h5py.File(filepath, "r")
tokenizer = {aa:i for i, aa in enumerate(list("ALGVSREDTIPKFQNYMHWC"))} #got from https://github.com/Rostlab/VESPA readme

def compute_model_logits(uniprot_id):
    dataset_name = "VESPAl"
    if uniprot_id in computed_data_file_handle.keys():
        data = computed_data_file_handle[uniprot_id][dataset_name] # h5py dataset object
        data = data[()] # numpy array: seq_len x 20, ALGVSREDTIPKFQNYMHWC got from https://github.com/Rostlab/VESPA readme
    else:
        data = np.array([]) # shape[0] is 0
    return data

def compute_variant_effect_scores(variants_df, prot_acc_version, output_logits):
    preds = []
    indices = variants_df[variants_df["prot_acc_version"]==prot_acc_version].index
    # print(len(indices))
    for idx in indices:
        tuple = variants_df.loc[idx]
        
        wt_tok_idx = tokenizer[tuple.wt]
        mt_tok_idx = tokenizer[tuple.mut]
        pos = tuple.prot_pos-1 #ncbi prot variants are 1 indexed, vespal preds are 0-indexed
        
        wt_logit = output_logits[pos][wt_tok_idx]
        mt_logit = output_logits[pos][mt_tok_idx]
        var_effect_score = mt_logit - wt_logit
        tuple = dict(tuple)
        tuple["pred"] = var_effect_score
        preds.append(tuple)
        # print(preds)
        # break
    return preds