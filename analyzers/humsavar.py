import sys
sys.path.append("../variant_effect_analysis")

from pyfaidx import Fasta
import pandas as pd


def extract_id(header):
    #print(header)
    return header.split('|')[1]

sequences = Fasta("data/uniprot_sprot.fasta", key_function=extract_id)
#print(sequences["P04637"])

f = open("data/humsavar.txt", "r")
for i, line in enumerate(f.readlines()):
    line = line.rstrip()

    gene_name = line[:10].strip()
    swissprot_ac = line[10:20].strip()
    ftid = line[20:31].strip()
    aa_change = line[31:45].strip()
    variant_category = line[45:53].strip()
    dbsnp = line[53:67].strip()
    disease_name = line[67:].strip()

    if swissprot_ac not in sequences:
        print(i, gene_name, swissprot_ac, "not found in uniprot")
        #break
    
    #print(i, gene_name, swissprot_ac, ftid, aa_change, variant_category, dbsnp, disease_name)
    #if i==30: break

