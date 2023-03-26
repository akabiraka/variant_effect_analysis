import sys
sys.path.append("../variant_effect_analysis")

import pandas as pd

def filter_drop_duplicates(df:pd.DataFrame, a_col_name:str):
    prev_num_rows = df.shape[0]
    # df[a_col_name].value_counts() # checking if some row occurs >twice
    # print("Number of rows that occurs >twice: ", df.groupby(a_col_name).filter(lambda x: len(x) > 1).shape[0])
    df = df.drop_duplicates(subset=a_col_name, keep="first")
    crnt_num_rows = df.shape[0]
    print(f"Number of duplicates rows removed: {prev_num_rows-crnt_num_rows}")
    return df


def filter_remove_null_nan_empty_entries(df:pd.DataFrame, a_col_name:str):
    prev_num_rows = df.shape[0]
    df = df[~pd.isna(df[a_col_name])]  # removing nan entries corresponding to a selected col column
    crnt_num_rows = df.shape[0]
    print(f"Number of NAN rows removed: {prev_num_rows-crnt_num_rows}")
    
    prev_num_rows = df.shape[0]    
    df = df[~pd.isnull(df[a_col_name])]  # removing null entries corresponding to a selected col column
    crnt_num_rows = df.shape[0]
    print(f"Number of NULL rows removed: {prev_num_rows-crnt_num_rows}")
    
    prev_num_rows = df.shape[0]
    df = df[df[a_col_name] != ""]
    crnt_num_rows = df.shape[0]
    print(f"Number of empty rows removed: {prev_num_rows-crnt_num_rows}")
    return df


def filter_multi_base_variants(x):
    # NC_000001.11:161362374:CTGG:ATGA
    return len(x)==1

def separate_ref_and_alt_allele(x):
    # NC_000023.11:154532268:C:A|NC_000023.11:154532045:A:C
    x = x.split("|")[0] 
    return x.split(":")[2], x.split(":")[3]