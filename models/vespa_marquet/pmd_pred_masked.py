import sys
sys.path.append("../variant_effect_analysis")
home_dir = ""

import time
import  pandas as pd

from models.aa_common.data_loader import get_pmd_dataset#, get_protein_sequences
import models.vespa_marquet.model_utils as model_utils


task = "pmd"
variants_df = get_pmd_dataset(home_dir)
variants_df = variants_df.rename(columns={"pmd_nr_id": "prot_acc_version"})
# protid_seq_tuple_list = get_protein_sequences(home_dir, max_seq_len=1022, return_type="protid_seq_tuple_list", data_type=task)

model_name = "vespa"
model, tokenizer = model_utils.get_model_tokenizer(model_name)
model_task_out_dir, model_logits_out_dir = model_utils.create_output_directories(model_name, task=task)

def execute(protid_seq_tuple_list):
    preds = []        
    for i, (prot_acc_version, mut_pos, seq) in enumerate(protid_seq_tuple_list):
        output_logits = model_utils.compute_model_logits_from_masked_sequences(model, tokenizer, prot_acc_version, seq, mut_pos, model_logits_out_dir)
        # print(output_logits.shape)
        preds_related_to_aprot = model_utils.compute_variant_effect_scores_from_masked_logits(variants_df, tokenizer, prot_acc_version, mut_pos, output_logits)
        # print(len(preds_related_to_aprot))
        preds += preds_related_to_aprot
    preds_df = pd.DataFrame(preds)   
    return preds_df

if __name__ == "__main__":
    start = time.time()
    prot_id_mutpos_seq_df = variants_df[["prot_acc_version", "prot_pos", "seq"]].drop_duplicates(keep="first") # prot_id_mutpos_seq_df
    data = list(prot_id_mutpos_seq_df.itertuples(index=False, name=None))

    chunk_size = 1 # 32 if torch.cuda.is_available() else 1
    data_chunks = [data[x:x+chunk_size] for x in range(0, len(data), chunk_size)]
    # data_chunks = data_chunks[:10] 
    print(f"#-of chunks: {len(data_chunks)}, 1st chunk size: {len(data_chunks[0])}")
    # print(data_chunks[0])
    # print(data_chunks[1])

    pred_dfs = []
    # sequential run and debugging
    # for i, data_chunk in enumerate(data_chunks):
    #     pred_df = execute(data_chunk)
    #     print(f"Finished {i}/{len(data_chunks)}th chunk: {pred_df.shape}")
    #     pred_dfs.append(pred_df)

    # mpi run    
    from mpi4py.futures import MPIPoolExecutor
    executor = MPIPoolExecutor()
    for i, pred_df in enumerate(executor.map(execute, data_chunks, unordered=True)):
        print(f"Finished {i}/{len(data_chunks)}th chunk: {pred_df.shape}")
        pred_dfs.append(pred_df)
    executor.shutdown()
    
    result_df = pd.concat(pred_dfs)  
    print("Saving predictions ...")
    result_df = result_df.rename(columns={"prot_acc_version": "pmd_nr_id"})
    result_df.to_csv(f"{model_task_out_dir}/preds_{model_name}_masked.tsv", sep="\t", index=False, header=True)
    print(result_df.shape)
    print(result_df.head())
        
    print(f"Time taken: {time.time()-start} seconds")