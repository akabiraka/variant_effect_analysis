import os
import sys
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)
home_dir = "../"

import time
import os
from Bio import Entrez
from urllib.error import HTTPError
import xml.etree.ElementTree as ET
import pandas as pd
import xmltodict


Entrez.email = "akabir0101@gmail.com" # provide your user email 
# RECOMMENDED: apply for API key from NCBI (https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/). 
# 10 queries per second with a valid API key, otherwise 3 queries per seconds are allowed for 'None'
Entrez.api_key = "328570309ccd040632796143ec88b51bcf08"
retmax = 500 # return 20 rs per batch example, max=1000
namespace = "https://www.ncbi.nlm.nih.gov/SNP/docsum"
ns = u'{%s}' % namespace
nsl = len(ns)

# dbSNP supported query terms (https://www.ncbi.nlm.nih.gov/snp/docs/entrez_help/) can be build and test online using web query builder (https://www.ncbi.nlm.nih.gov/snp/advanced) 
# esearch handle
eShandle = Entrez.esearch(db="snp",  # search dbSNP
                          term='"humo sapiens"[Organism] AND "missense variant"[Function Class] AND "snp protein"[Filter] AND "validated by alfa"[Filter] ',
                          usehistory="y", #cache result on server for download in batches
                          retmax=retmax
                         )                        

# get esearch result
eSresult = Entrez.read(eShandle)
webenv = eSresult["WebEnv"]
query_key = eSresult["QueryKey"]
total_count = int(eSresult["Count"])
print(f"Query result count:: {total_count}, Fetch count: {len(range(0, total_count, retmax))}")


def parse_response_xml(x:str):
    o = xmltodict.parse(x) #returns dict obj
    docs = o["ExchangeSet"]["DocumentSummary"]

    outs = []
    for i, doc in enumerate(docs):
        genes = doc["GENES"]["GENE_E"] # can have multiples
        if isinstance(genes, list):
            genes = [f'{gene["NAME"]}:{gene["GENE_ID"]}' for gene in genes]
            genes = ",".join(genes)
        elif isinstance(genes, dict):
            genes = f'{genes["NAME"]}:{genes["GENE_ID"]}'
        # print(genes)

        
        mafs = doc["GLOBAL_MAFS"]["MAF"] # can have multiples
        if isinstance(mafs, list):
            mafs = [f'{maf["STUDY"]}:{maf["FREQ"]}' for maf in mafs]
            mafs = ",".join(mafs)
        elif isinstance(mafs, dict):
            mafs = f'{mafs["STUDY"]}:{mafs["FREQ"]}'
            # print(mafs)

        variations = doc["DOCSUM"] # an example: HGVS=NC_000023.11:g.154360389T>A,NC_000023.11:g.154360389T>C,NW_003871103.3:g.1794368T>A,NW_003871103.3:g.1794368T>C,NG_011506.2:g.19250A>T,NG_011506.2:g.19250A>G,NM_001456.4:c.3406A>T,NM_001456.4:c.3406A>G,NM_001456.3:c.3406A>T,NM_001456.3:c.3406A>G,NM_001110556.2:c.3406A>T,NM_001110556.2:c.3406A>G,NM_001110556.1:c.3406A>T,NM_001110556.1:c.3406A>G,NC_000023.10:g.153588757T>A,NC_000023.10:g.153588757T>C,NP_001447.2:p.Ile1136Phe,NP_001447.2:p.Ile1136Val,NP_001104026.1:p.Ile1136Phe,NP_001104026.1:p.Ile1136Val|SEQ=[T/A/C]|LEN=1|GENE=FLNA:2316
        variations = (variations[5:].split("|")[0]).split(",")
        variations = [v for v in variations if v.startswith("NP_")]
        variations = ",".join(variations)
        
        data = {
            "snp_id": doc["SNP_ID"],
            "acc": doc["ACC"],
            "chrpos": doc["CHRPOS"],
            "spdi": doc["SPDI"],
            "tax_id": doc["TAX_ID"],
            "snp_class": doc["SNP_CLASS"],
            "create_date": doc["CREATEDATE"],
            "update_date": doc["UPDATEDATE"],
            "clinical_significance": doc["CLINICAL_SIGNIFICANCE"],
            "fxn_class": doc["FXN_CLASS"],
            "validated": doc["VALIDATED"],
            "genes": genes,
            "mafs": mafs,
            "variations": variations
        }
        outs.append(data)
        
        # if i==5: break
    return outs

def save(data, out_filepath:str):
    # data: list of objects
    out_df = pd.DataFrame(data)
    if not os.path.exists(out_filepath):
        out_df.to_csv(out_filepath, chunksize=10000, sep="\t", index=False, mode="a", header=True)
    else:
        out_df.to_csv(out_filepath, chunksize=10000, sep="\t", index=False, mode="a", header=False)
        
        

def download_batch(start, retmax):
    attempt = 0
    while (attempt < 3):
        attempt += 1
        try:
            fetch_handle = Entrez.efetch(db="snp",
                                        # rettype="",
                                        retmode="xml",
                                        retstart=start,
                                        retmax=retmax,
                                        webenv=webenv,
                                        query_key=query_key)
            break
        except HTTPError as err:
            if 400 <= err.code <= 599:
                print("Received error from server %s" % err)
                print("Attempt %i of 3" % attempt)
                time.sleep(10)
            else:
                raise
    try:                
        data = fetch_handle.read().decode()
        data = parse_response_xml(data)
        fetch_handle.close()
        return data
    except:
        print(f"Error. Downloading again record {start+1} to {start+retmax}")
        return download_batch(start, retmax)
    # return fetch_handle
    
    
    
# sample codes adopted with modifications from http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc139.
out_filepath = "data/dbsnp/search_results/dbsnp_HumanMissenseALFAVariants_3.txt"
# if os.path.exists(out_filepath): os.remove(out_filepath)


fetch_count = 8000 # 1-indexed
start_idx = 4000000 # 0-indexed
for start in range(start_idx, total_count, retmax):
    end = min(total_count, start+retmax)
    print(f"Fetch no: {fetch_count}\tDownloading record {start+1} to {end}")
    
    # fetch_handle = download_batch(start, retmax)
    # if (fetch_handle):
    #     data = fetch_handle.read().decode()
    #     data = parse_response_xml(data)
    #     # print(data)
    #     save(data, out_filepath)
    #     fetch_handle.close()
        
    data = download_batch(start, retmax)
    save(data, out_filepath)
        
    # if fetch_count==2: break
    fetch_count += 1            