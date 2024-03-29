import sys
sys.path.append("../variant_effect_analysis")
home_dir = ""

import os
import torch
import numpy as np
import pandas as pd
from models.aa_common.data_loader import get_pmd_dbnsfp_dataset
import models.tape_rao.model_utils as model_utils
import utils.pickle_utils as pickle_utils
from tape import ProteinBertModel, TAPETokenizer

task = "pmd"
variants_df, protid_seq_dict = get_pmd_dbnsfp_dataset(home_dir)

model_name = "protbert"
model = ProteinBertModel.from_pretrained('bert-base')
tokenizer = TAPETokenizer(vocab='iupac') 
model_task_out_dir, model_logits_out_dir = model_utils.create_output_directories(model_name, task, home_dir)


def get_embedding(seq, filename):
    filepath = f"{model_logits_out_dir}{filename}.pkl"

    if os.path.exists(filepath):
        # print(f"Model logits already exists: {filename}")
        embedding = pickle_utils.load_pickle(filepath) 
    else: 
        # print(f"Computing model logits: {filename}")
        with torch.no_grad():
            token_ids = torch.tensor(np.array([tokenizer.encode(seq)]))
            embedding = model(token_ids)[0].squeeze(0).detach().numpy()
        pickle_utils.save_as_pickle(embedding, filepath)
    return embedding


def compute_variant_effect_score(protid, seq, one_indexed_mut_pos, wt_aa, mt_aa):
    wt_seq = list(seq)
    mt_seq = list(seq)
    mt_seq[one_indexed_mut_pos-1] = mt_aa

    wt_filename = f"{protid}"
    mt_filename = f"{protid}_{str(one_indexed_mut_pos)}_{mt_aa}"
    
    wt_embedding = get_embedding(wt_seq, wt_filename)[1:-1] # 1st and last tokens are <cls>=2 and <sep>=3
    mt_embedding = get_embedding(mt_seq, mt_filename)[1:-1]
    # print("Computing mut-type seq embedding: ", mt_filename)
    # print(wt_embedding.shape, mt_embedding.shape)

    effect_score = abs(mt_embedding - wt_embedding).sum() / (768*len(seq)) # embedding_dim = 768
    # print(effect_score)
    return effect_score


def execute(row):
    index = row[0]
    data = row[1]
    
    protid, seq, one_indexed_mut_pos, wt_aa, mt_aa = data["prot_acc_version"], protid_seq_dict[data["prot_acc_version"]], data["prot_pos"], data["wt"], data["mut"]
    effect_score = compute_variant_effect_score(protid, seq, one_indexed_mut_pos, wt_aa, mt_aa)
    print(index, protid, one_indexed_mut_pos, wt_aa, mt_aa)

    row = variants_df.loc[index]
    row = dict(row)
    row["pred"] = effect_score
    
    return row

def compute_wt_embedding(i, protid_seq_tuple):
    protid, seq = protid_seq_tuple
    get_embedding(seq, f"{protid}")
    print(f"Computing wild-type seq embedding: {i}|{len(protid_seq_dict)}: {protid}")

if __name__ == "__main__":

    # variants_df = variants_df.head(10)
    preds = [] 
    
    from mpi4py.futures import MPIPoolExecutor
    executor = MPIPoolExecutor()
    executor.map(compute_wt_embedding, list(range(len(protid_seq_dict))), list(protid_seq_dict.items()), unordered=False)
    executor.shutdown()

    executor = MPIPoolExecutor()
    for i, row in enumerate(executor.map(execute, variants_df.iterrows(), unordered=False)):
        preds.append(row)
    executor.shutdown()

    result_df = pd.DataFrame(preds)   
    print(result_df.shape)
    print(result_df.head())
    result_df.to_csv(f"{model_task_out_dir}/preds_{model_name}_embed.tsv", sep="\t", index=False, header=True)