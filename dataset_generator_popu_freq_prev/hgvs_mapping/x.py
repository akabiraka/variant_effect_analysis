import os
import sys
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)
home_dir = ""

import time
import os
import pandas as pd

import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.assemblymapper
import hgvs.normalizer


hp = hgvs.parser.Parser()
# initialize the mapper for GRCh38 with splign-based alignments
hdp = hgvs.dataproviders.uta.connect()
am38 = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name='GRCh38', alt_aln_method='splign', replace_reference=True)
hn = hgvs.normalizer.Normalizer(hdp)


def map_genomic_variants_to_coding_variants(row):
    # print(row)
    # raise
    hgvs_transcript_prot_variant_pairs = []
    for alt_allele in row.alfa_alt_alleles.split(','):
        alt_allele = alt_allele.strip()

        hgvs_g = f"{row.alfa_chrom_acc_version}:g.{row.alfa_1indexed_chrom_pos}{row.alfa_ref_allele}>{alt_allele}" # ie "NW_003571061.2:g.275107T>C" g. for genomic variants
        var_g = hn.normalize(hp.parse_hgvs_variant(hgvs_g))

        transcripts = am38.relevant_transcripts(var_g)
        mane_coding_transcript = row.mane_refseq_nuc  

        hgvs_coding_transcript = None
        if mane_coding_transcript in transcripts:
            hgvs_coding_transcript = mane_coding_transcript
        else:
            for t in transcripts:
                if t.split(".")[0] == mane_coding_transcript.split(".")[0]: # matching transcript acc only
                    hgvs_coding_transcript = t
                    break

        print("\t", row.snp_id, hgvs_coding_transcript, str(var_g), end=" ")

        if hgvs_coding_transcript is not None:
            var_c = am38.g_to_c(var_g, hgvs_coding_transcript)
            var_p = am38.c_to_p(var_c)
            var_p.posedit.uncertain = False
            print(str(var_c), str(var_p))
            
            hgvs_transcript_prot_variant_pairs.append(hgvs_coding_transcript+"-"+str(var_p))
    return ",".join(hgvs_transcript_prot_variant_pairs)
    # raise


def run_df(df:pd.DataFrame, chunk_id=0):
    start_time = time.time()
    print("starting:", chunk_id)
    out_filepath = home_dir+f"out_dir/snps_with_mane_alfa_protvariants_{chunk_id}.tsv"
    if not os.path.exists(out_filepath):
        df["hgvs_transcript_prot_variant_pairs"] = df.apply(map_genomic_variants_to_coding_variants, axis=1)
        df.to_csv(out_filepath, sep="\t", index=False, header=True)
    print("finished:", chunk_id, time.time()-start_time, "seconds to exe")


# inp_filepath = home_dir+"snps_with_mane_alfa_short_indexed.tsv"
# snps_with_mane_alfa_df_iterator = pd.read_csv(inp_filepath, chunksize=10, sep="\t", index_col=0)

inp_filepath = home_dir+"snps_with_mane_alfa_short.tsv"
snps_with_mane_alfa_df_iterator = pd.read_csv(inp_filepath, chunksize=10, sep="\t")

for i, snps_with_mane_alfa_chunk_df in enumerate(snps_with_mane_alfa_df_iterator):
    print(snps_with_mane_alfa_chunk_df)
    run_df(snps_with_mane_alfa_chunk_df, i)
    if i==1: break


# from joblib import Parallel, delayed
# Parallel(n_jobs=2)(delayed(run_df)(snps_with_mane_alfa_chunk_df, i) 
#                    for i, snps_with_mane_alfa_chunk_df in enumerate(snps_with_mane_alfa_df_iterator) if i<2)


# x = [i for i, snps_with_mane_alfa_chunk_df in enumerate(snps_with_mane_alfa_df_iterator) if i<2]
# print(x)