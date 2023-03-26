import sys
sys.path.append("../variant_effect_analysis")

import os
import time
import gzip
import pandas as pd
import multiprocessing as mp


def get_column_names(vcf_path):
    with gzip.open(vcf_path, "rt") as ifile:
        for line in ifile:
            #print(line)
            if line.startswith("#CHROM"):
                line = line[1:]
                line = line.rstrip()
                vcf_names = [x for x in line.split('\t')]
                break
    ifile.close()
    return vcf_names


def print_SNP_summary(vcf_chunk_iterator):
    for chunk_df in vcf_chunk_iterator:
        print(chunk_df.shape)
        # print(chunk_df.head())
        # print(chunk_df.loc[0])
        break


def is_valid_alt_allele_count(alt_allele_count):
    if alt_allele_count > 0: # alt_allele_count must be >0
        return True
    else: return False



def run_each_chunk(chunk_df:pd.DataFrame):
    start = time.time()
    out_df = pd.DataFrame(columns=chunk_df.columns)
    for row_no, row in enumerate(chunk_df.itertuples()):
        # print(row)
        alts = row.ALT.split(",")
        total_allele_count, alt_allele_counts =  row.SAMN10492705.split(":")
        total_allele_count, alt_allele_counts = int(total_allele_count), list(map(int, alt_allele_counts.split(",")))
        total_allele_count = 1 if total_allele_count==0 else total_allele_count
        allele_freqs = [x / total_allele_count for x in alt_allele_counts]
        #print(allele_freqs)

        for i, alt in enumerate(alts):
            # add more filters here
            if is_valid_alt_allele_count(alt_allele_counts[i]):
                new_row = row._replace(SAMN10492705=str(total_allele_count)+":"+str(alt_allele_counts[i]))
                new_row = new_row._replace(ALT=alt)
                new_row = pd.Series(new_row._asdict()).to_frame().T
                out_df = pd.concat([out_df, new_row], ignore_index=True)
                # print(new_row)
        
        # if is_valid_allele_freq_filter(allele_freqs): 
        #     new_row = pd.Series(row._asdict()).to_frame().T
        #     out_df = pd.concat([out_df, new_row], ignore_index=True)
        #     # print(chunk_no, row_no, new_row, out_df)

        
        # if row_no==20: break
    out_df.to_csv(out_filepath, compression='gzip', chunksize=10000, sep="\t", index=False, mode="a")
    end = time.time()
    exe_time = (end-start)*10**3 #in miliseconds 
    return exe_time

filepath = "data/snp_population_freq/2020_11.vcf.gz"
column_names = get_column_names(filepath)
data_chunk_iterator = pd.read_csv(filepath, compression='gzip', comment='#', chunksize=10000, delim_whitespace=False, sep="\t", header=None, names=column_names)
#print_SNP_summary(data_chunk_iterator)

out_filepath = "data/snp_population_freq/2020_11_filtered.vcf.gz"
if os.path.exists(out_filepath): os.remove(out_filepath)



# n_chunks_to_skip = 0
# for chunk_no, chunk_df in enumerate(data_chunk_iterator):
#     if chunk_no < n_chunks_to_skip: continue
#     start = time.time()
#     # out_df = run_each_chunk(chunk_no, chunk_df)
      # out_df.to_csv(out_filepath, compression='gzip', chunksize=10000, sep="\t", index=False, mode="a")  
#     end = time.time()
#     print(f"Chunk: {chunk_no}\tTime taken: {(end-start)*10**3:.03f}ms\t{out_df.shape}")
    
#     if chunk_no==2: break

n_workers = mp.cpu_count() 
print(n_workers)
processes = mp.Pool(n_workers)
tasks = processes.imap(run_each_chunk, data_chunk_iterator, chunksize=100)
for task_out in tasks:
    exe_time = task_out
    print(f"Chunk: \tTime taken: {exe_time:.03f}ms")
    
    # test filtered dataframe
    # filtered_chunk_iterator = pd.read_csv(out_filepath, compression='gzip', comment='#', chunksize=10000, delim_whitespace=False, sep="\t") # "Index" col stands for the indices of the original data
    # print_SNP_summary(filtered_chunk_iterator)
