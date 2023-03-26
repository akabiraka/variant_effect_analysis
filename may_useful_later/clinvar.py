import sys
sys.path.append("../variant_effect_analysis")

import pandas as pd


def print_variant_summary_data(filepath):
    df = pd.read_csv(filepath, sep="\t")

    print(df.shape)
    print(df.head())
    print(df.columns)
    print(df.loc[0])
    
    print(len(df[(df["ClinSigSimple"]==0) & (df["Assembly"]=="GRCh38")])) # #-of rows where ClinSigSimple=0 and Assembly=GRCh38
    print(len(df[(df["ClinSigSimple"]==1) & (df["Assembly"]=="GRCh38")])) # #-of rows where ClinSigSimple=1 and Assembly=GRCh38
    print(len(df[(df["ClinSigSimple"]==-1) & (df["Assembly"]=="GRCh38")])) # #-of rows where ClinSigSimple=-1 and Assembly=GRCh38

#filepath = "data/clinvar/2021_12/variant_summary.txt"
filepath = "data/clinvar/2023_01/variant_summary.txt"
print_variant_summary_data(filepath)
