import sys
home_dir = ""
sys.path.append("../variant_effect_analysis")

import time
import pandas as pd

# from models.aa_common.data_loader import get_protein_sequences, get_pathogenicity_analysis_SNVs
from models.aa_common.data_loader import get_patho_and_likelypatho_SNVs, generate_neutral_SNVs, get_protein_sequences
import models.tape_rao.model_utils as model_utils

list_of_neutral_variants_for_patho_df = generate_neutral_SNVs(home_dir, pathogenicity_type="pathogenic")
list_of_neutral_variants_for_likelypatho_df = generate_neutral_SNVs(home_dir, pathogenicity_type="likely_pathogenic")

patho_and_likelypatho_variants_df = get_patho_and_likelypatho_SNVs(home_dir)
protid_seq_tuple_list = get_protein_sequences(home_dir=home_dir, max_seq_len=1022, return_type="protid_seq_tuple_list", data_type="patho_and_likelypatho")

model_name = "unirep" # unirep, protbert
model, tokenizer = model_utils.get_model_tokenizer(model_name)
model_task_out_dir, model_logits_out_dir, patho_out_dir, likelypatho_out_dir = model_utils.create_output_directories(model_name)

# def execute(protid_seq_tuple_list, analysis_no=0):
#     variants_df = list_of_variants_df[analysis_no]
#     preds = []        
#     for i, (prot_acc_version, seq) in enumerate(protid_seq_tuple_list):
#         output_logits = model_utils.compute_model_logits(model, tokenizer, prot_acc_version, seq, model_logits_out_dir)
#         preds_related_to_aprot = model_utils.compute_variant_effect_scores(variants_df, tokenizer, prot_acc_version, output_logits)
#         preds += preds_related_to_aprot
#     preds_df = pd.DataFrame(preds)   
#     return preds_df

def execute_neutral_df(variants_df, prot_acc_version, output_logits):
    preds_related_to_aprot = model_utils.compute_variant_effect_scores(variants_df, tokenizer, prot_acc_version, output_logits)
    return preds_related_to_aprot


def execute(protid_seq_tuple_list):
    patho_and_likelypatho_preds = []
    neutral_patho_preds = [[] for i in range(10)]
    neutral_likelypatho_preds = [[] for i in range(10)]
    for i, (prot_acc_version, seq) in enumerate(protid_seq_tuple_list):
        output_logits = model_utils.compute_model_logits(model, tokenizer, prot_acc_version, seq, model_logits_out_dir)

        preds_related_to_aprot = model_utils.compute_variant_effect_scores(patho_and_likelypatho_variants_df, tokenizer, prot_acc_version, output_logits)
        patho_and_likelypatho_preds += preds_related_to_aprot

        for analysis_no in range(10):
            neutral_patho_preds[analysis_no] += execute_neutral_df(list_of_neutral_variants_for_patho_df[analysis_no], prot_acc_version, output_logits)
            neutral_likelypatho_preds[analysis_no] += execute_neutral_df(list_of_neutral_variants_for_likelypatho_df[analysis_no], prot_acc_version, output_logits)

    patho_and_likelypatho_preds_df = pd.DataFrame(patho_and_likelypatho_preds)   
    neutral_patho_preds_dfs = [pd.DataFrame(pred) for pred in neutral_patho_preds]
    neutral_likelypatho_preds_dfs = [pd.DataFrame(pred) for pred in neutral_likelypatho_preds]

    return patho_and_likelypatho_preds_df, neutral_patho_preds_dfs, neutral_likelypatho_preds_dfs


if __name__ == "__main__":
    start = time.time()
    
    data = protid_seq_tuple_list

    chunk_size = 1 # 32 if torch.cuda.is_available() else 1
    data_chunks = [data[x:x+chunk_size] for x in range(0, len(data), chunk_size)]
    # data_chunks = data_chunks[:10] 
    print(f"#-of chunks: {len(data_chunks)}, 1st chunk size: {len(data_chunks[0])}")
    
    patho_and_likelypatho_pred_dfs = []
    list_of_neutral_patho_preds_dfs = [[] for i in range(10)]
    list_of_neutral_likelypatho_preds_dfs = [[] for i in range(10)]

    # sequential run and debugging
    for i, data_chunk in enumerate(data_chunks):
        patho_and_likelypatho_preds_df, neutral_patho_preds_dfs, neutral_likelypatho_preds_dfs = execute(data_chunk)
        
        print(f"Finished {i}/{len(data_chunks)}th chunk: {patho_and_likelypatho_preds_df.shape}")
        patho_and_likelypatho_pred_dfs.append(patho_and_likelypatho_preds_df)
        for i in range(10):
            list_of_neutral_patho_preds_dfs[i].append(neutral_patho_preds_dfs[i]) 
            list_of_neutral_likelypatho_preds_dfs[i].append(neutral_likelypatho_preds_dfs[i]) 

    result_df = pd.concat(patho_and_likelypatho_pred_dfs)  
    print("Saving predictions ...")
    result_df.to_csv(f"{model_task_out_dir}/patho_and_likelypatho.csv", sep="\t", index=False, header=True)
    print(result_df.shape)
    # print(result_df.head())

    for i in range(10):
        result_df = pd.concat(list_of_neutral_patho_preds_dfs[i])  
        print("Saving predictions ...")
        result_df.to_csv(f"{patho_out_dir}/{i}.csv", sep="\t", index=False, header=True)
        print(result_df.shape)

    for i in range(10):
        result_df = pd.concat(list_of_neutral_likelypatho_preds_dfs[i])  
        print("Saving predictions ...")
        result_df.to_csv(f"{likelypatho_out_dir}/{i}.csv", sep="\t", index=False, header=True)
        print(result_df.shape)

        
    # for analysis_no in range(10):
    #     pred_dfs = []
    #     # sequential run and debugging
    #     for i, data_chunk in enumerate(data_chunks):
    #         pred_df = execute(data_chunk, analysis_no)
    #         print(f"Finished {i}/{len(data_chunks)}th chunk: {pred_df.shape}")
            # pred_dfs.append(pred_df)

        # mpi run    
        # from mpi4py.futures import MPIPoolExecutor
        # executor = MPIPoolExecutor()
        # for i, pred_df in enumerate(executor.map(execute, data_chunks, unordered=True)):
        #     # print(f"Finished {i}/{len(data_chunks)}th chunk: {pred_df.shape}")
        #     pred_dfs.append(pred_df)
        # executor.shutdown()
        
    #     result_df = pd.concat(pred_dfs)  
    #     print("Saving predictions ...")
    #     result_df.to_csv(f"{model_task_out_dir}/{str(analysis_no)}.csv", sep="\t", index=False, header=True)
    #     print(result_df.shape)
    #     # print(result_df.head())
    #     break
        
    # print(f"Time taken: {time.time()-start} seconds")