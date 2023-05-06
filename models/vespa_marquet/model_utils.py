import sys
sys.path.append("../variant_effect_analysis")

import re
import os
import numpy as np
import torch

import utils.pickle_utils as pickle_utils
import h5py

from vespa.predict.logodds import T5_condProbas

def get_model_tokenizer(home_dir=""):
    print("\nLog: Model loading ...")
    cache_dir="models/vespa_marquet/cache"
    t5_condProbas = T5_condProbas(cache_dir=cache_dir)
    model, tokenizer = t5_condProbas.prott5.get_model(1) #1=LOGODDS
    return model, tokenizer

def compute_model_logits(model, tokenizer, prot_acc_version, seq, logits_output_path)->np.array:
    filepath = f"{logits_output_path}{prot_acc_version}.pkl"
    if os.path.exists(filepath):
        print(f"Model logits already exists: {prot_acc_version}")
        logits = pickle_utils.load_pickle(filepath) 
    else: 
        print(f"Computing model logits: {prot_acc_version}")
        with torch.no_grad():
            seq = re.sub(r"[UZOB]", "X", seq) # replacing unknown amino acids 
            seq = " ".join(list(seq))
            input_ids = tokenizer(seq, return_tensors="pt").input_ids.to("cpu")
            logits = model(input_ids, labels=input_ids).logits
            logits = logits[0].detach().numpy() 
            pickle_utils.save_as_pickle(logits, filepath)
    # print(logits.shape)
    return logits # logits shape: (seq_len+1, vocab_size=128) # eos token added

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
static_tokenizer = {aa:i for i, aa in enumerate(list("ALGVSREDTIPKFQNYMHWC"))} #got from https://github.com/Rostlab/VESPA readme

def get_model_logits(uniprot_id):
    dataset_name = "VESPAl"
    if uniprot_id in computed_data_file_handle.keys():
        data = computed_data_file_handle[uniprot_id][dataset_name] # h5py dataset object
        data = data[()] # numpy array: seq_len x 20, ALGVSREDTIPKFQNYMHWC got from https://github.com/Rostlab/VESPA readme
    else:
        data = np.array([]) # shape[0] is 0
    return data

def compute_variant_effect_scores_vespal(variants_df, prot_acc_version, output_logits):
    """This is for population freq and pathogenicity analysis
    """
    preds = []
    indices = variants_df[variants_df["prot_acc_version"]==prot_acc_version].index
    # print(len(indices))
    for idx in indices:
        tuple = variants_df.loc[idx]
        
        wt_tok_idx = static_tokenizer[tuple.wt]
        mt_tok_idx = static_tokenizer[tuple.mut]
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


def compute_variant_effect_scores(variants_df, tokenizer, prot_acc_version, output_logits, model_aa_prefix="▁"):
    """this is robust when computing from the model directly
    """
    preds = []
    indices = variants_df[variants_df["prot_acc_version"]==prot_acc_version].index 
    for idx in indices:
        tuple = variants_df.loc[idx]
        
        wt_tok_idx = tokenizer.convert_tokens_to_ids(model_aa_prefix+tuple.wt)
        mt_tok_idx = tokenizer.convert_tokens_to_ids(model_aa_prefix+tuple.mut)
        pos = tuple.prot_pos-1 #ncbi prot variants are 1 indexed, outputs are 0-indexed, so -1

        # print(output_logits.shape, pos)
        if pos>=output_logits.shape[0]: continue
        
        wt_logit = output_logits[pos][wt_tok_idx]
        mt_logit = output_logits[pos][mt_tok_idx]
        var_effect_score = mt_logit - wt_logit
        tuple = dict(tuple)
        tuple["pred"] = var_effect_score
        preds.append(tuple)
        # print(preds)
        # break
    return preds