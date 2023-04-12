import sys
home_dir = ""
sys.path.append("../variant_effect_analysis")

import torch
import pandas as pd
import time

from models.aa_common.data_loader import get_population_freq_SNVs, get_protein_sequences
import models.bioembeddings_dallago.model_utils as model_utils

task = "popu_freq"
variants_df = get_population_freq_SNVs()
protid_seq_tuple_list = get_protein_sequences(home_dir=home_dir, max_seq_len=1022, return_type="protid_seq_tuple_list", data_type=task)

model_name = "prottrans_t5_bfd" #  plus_rnn, prottrans_bert_bfd, prottrans_albert_bfd, prottrans_xlnet_uniref100, prottrans_t5_bfd, prottrans_t5_uniref50, prottrans_t5_xl_u50
model, tokenizer, model_name = model_utils.get_model_tokenizer(model_name) 
model_task_out_dir, model_logits_out_dir = model_utils.create_output_directories(model_name, task=task)

def execute(protid_seq_tuple_list):
    preds = []        
    for i, (prot_acc_version, seq) in enumerate(protid_seq_tuple_list):
        output_logits = model_utils.compute_model_logits(model, prot_acc_version, seq, model_logits_out_dir)
        preds_related_to_aprot = model_utils.compute_variant_effect_scores(variants_df, tokenizer, prot_acc_version, output_logits, model.aa_prefix)
        preds += preds_related_to_aprot
    preds_df = pd.DataFrame(preds)   
    return preds_df


if __name__=="__main__": # main worker
    start = time.time()

    data = protid_seq_tuple_list

    chunk_size = 1 # 32 if torch.cuda.is_available() else 1
    data_chunks = [data[x:x+chunk_size] for x in range(0, len(data), chunk_size)]
    # data_chunks = data_chunks[:10] 
    print(f"#-of chunks: {len(data_chunks)}, 1st chunk size: {len(data_chunks[0])}")

    pred_dfs = []
    # # sequential run and debugging
    for i, data_chunk in enumerate(data_chunks):
        pred_df = execute(data_chunk)
        print(f"Finished {i}/{len(data_chunks)}th chunk: {pred_df.shape}")
        pred_dfs.append(pred_df)

    # mpi run    
    # from mpi4py.futures import MPIPoolExecutor
    # executor = MPIPoolExecutor()
    # for i, pred_df in enumerate(executor.map(execute, data_chunks, unordered=True)):
    #     print(f"Finished {i}/{len(data_chunks)}th chunk: {pred_df.shape}")
    #     pred_dfs.append(pred_df)
    # executor.shutdown()

    result_df = pd.concat(pred_dfs)  
    print("Saving predictions ...")  
    result_df.to_csv(f"{model_task_out_dir}/preds_{model_name}.csv", sep="\t", index=False, header=True)
    print(result_df.shape)
    print(result_df.head())

    print(f"Time taken: {time.time()-start} seconds")