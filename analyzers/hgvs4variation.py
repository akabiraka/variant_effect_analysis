import sys
sys.path.append("../variant_effect_analysis")

import gzip
import pandas as pd

def get_col_names(filepath):
    with gzip.open(filepath, "rt") as f:
        for line_no, line in enumerate(f):
            # print(line)
            if line.startswith("#Symbol"):
                line = line[1:] # skipping the # symbol
                line = line.rstrip()
                col_names = [x for x in line.split('\t')]
                # if line_no==100: break
                break
        return col_names

def print_summary(dfs_iterator:iter):
    first_df = dfs_iterator.__next__()
    print(first_df.shape)
    print(first_df.head())
    print(first_df.loc[0])
        
filepath = "data/clinvar/2023_01/hgvs4variation.txt.gz"
col_names = get_col_names(filepath)
# hgvs_dfs_iterator = pd.read_csv(filepath, compression='gzip', comment='#', chunksize=10000, delim_whitespace=False, sep="\t", header=None, names=col_names)
# # print_summary(dfs_iterator)

# for df_idx, hgvs_df in enumerate(hgvs_dfs_iterator):
#     hits = hgvs_df[hgvs_df["VariationID"]==3]
    
#     if hits.shape[0]>0: 
#         print(df_idx, hits)
#         break

hgvs_df = pd.read_csv(filepath, compression='gzip', comment='#', delim_whitespace=False, sep="\t", header=None, names=col_names)    
print(hgvs_df.shape)

print(hgvs_df[hgvs_df["VariationID"]==1320032])