import sys
home_dir = ""
sys.path.append("../variant_effect_analysis")

import time
import pandas as pd

from models.aa_common.data_loader import get_patho_and_likelypatho_SNVs, get_protein_sequences
import models.vespa_marquet.model_utils as model_utils

task = "patho_and_likelypatho"
patho_and_likelypatho_variants_df = get_patho_and_likelypatho_SNVs(home_dir)

protid_seq_tuple_list = get_protein_sequences(home_dir=home_dir, max_seq_len=1022, return_type="protid_seq_tuple_list", data_type=task)

model_name = "vespal"
model_task_out_dir, model_logits_out_dir = model_utils.create_output_directories(model_name, task, home_dir)

# np_to_uniprot_mapping_df = pd.read_csv(home_dir+"data/gene/np_to_uniprot_mapping.csv", sep="\t")
uniprot_to_np_mapping_df = pd.read_csv(home_dir+"models/vespa_marquet/cache/uniprot_to_refseq_mapping.tsv", sep="\t")

def execute(protid_seq_tuple_list):
    variants_df = patho_and_likelypatho_variants_df.copy(deep=True)
    preds = []   
    for i, (prot_acc_version, seq) in enumerate(protid_seq_tuple_list):

        preds_related_to_aprot = None
        np_uniprot_pairs_df = uniprot_to_np_mapping_df[uniprot_to_np_mapping_df["refseq_id"]==prot_acc_version].reset_index()
        print(prot_acc_version, np_uniprot_pairs_df)

        if np_uniprot_pairs_df.shape[0]== 1: # if NP-id mapped to 1 uniprot id
            uniprot_id = np_uniprot_pairs_df.at[0, "uniprot_id"]
            # print(uniprot_id)
            output_logits = model_utils.get_model_logits(uniprot_id)
            if output_logits.shape[0] == len(seq): # corresponding to same sequence length
                preds_related_to_aprot = model_utils.compute_variant_effect_scores(variants_df, prot_acc_version, output_logits)
        
        else:
            for tuple in np_uniprot_pairs_df.itertuples():
                uniprot_id = tuple.uniprot_id
                output_logits = model_utils.get_model_logits(uniprot_id)
                
                # print(output_logits.shape[0], len(seq))
                if output_logits.shape[0] == len(seq):
                    preds_related_to_aprot = model_utils.compute_variant_effect_scores(variants_df, prot_acc_version, output_logits)
                    # print(uniprot_id)
                    break
        if preds_related_to_aprot is not None: 
            preds += preds_related_to_aprot

    preds_df = pd.DataFrame(preds)   
    return preds_df
        # break


if __name__ == "__main__":
    start = time.time()
    
    data = protid_seq_tuple_list

    chunk_size = 1 # 32 if torch.cuda.is_available() else 1
    data_chunks = [data[x:x+chunk_size] for x in range(0, len(data), chunk_size)]
    # data_chunks = data_chunks[:10] 
    print(f"#-of chunks: {len(data_chunks)}, 1st chunk size: {len(data_chunks[0])}")
        
    pred_dfs = []
    # sequential run and debugging
    for i, data_chunk in enumerate(data_chunks):
        # if i<475: continue
        pred_df = execute(data_chunk)
        print(f"Finished {i}/{len(data_chunks)}th chunk: {pred_df.shape}")
        pred_dfs.append(pred_df)
        # if i==20: break

        # mpi run    
        # from mpi4py.futures import MPIPoolExecutor
        # executor = MPIPoolExecutor()
        # for i, pred_df in enumerate(executor.map(execute, data_chunks, unordered=True)):
        #     # print(f"Finished {i}/{len(data_chunks)}th chunk: {pred_df.shape}")
        #     pred_dfs.append(pred_df)
        # executor.shutdown()
        
    result_df = pd.concat(pred_dfs)  
    print("Saving predictions ...")
    result_df.to_csv(f"{model_task_out_dir}/preds_{model_name}.csv", sep="\t", index=False, header=True)
    print(result_df.shape)
    # print(result_df.head())
        
    print(f"Time taken: {time.time()-start} seconds")