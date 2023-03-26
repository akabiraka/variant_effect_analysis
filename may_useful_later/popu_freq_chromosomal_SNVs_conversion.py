import sys
sys.path.append("../variant_effect_analysis")
home_dir = ""

import pandas as pd


def get_chromosomal_variants_df(dbsnps_df):
    chromosomal_variants = []
    for i, tuple in enumerate(dbsnps_df.itertuples()):

        if len(tuple.REF)>1: # only considering single neucleodite variants
            continue
        
        wt_population, mut_poulations = int(tuple.SAMN10492705.split(":")[0]), tuple.SAMN10492705.split(":")[1].split(",")
        total_population = wt_population+sum(list(map(int, mut_poulations)))


        alt_alleles = tuple.ALT.split(",")
        try:
            for j, alt_allele in enumerate(alt_alleles):
                if len(alt_allele) > 1: continue # only considering SNVs (single nucleotide variants)

                if j < len(mut_poulations):
                    mut_poulation = int(mut_poulations[j])
                else: mut_poulation = 0

                v = {"snp_id": tuple.snp_id,
                    "chrom_acc_version": tuple.CHROM,
                    "pos": tuple.POS, # 1 indexed
                    "ref": tuple.REF,
                    "alt": alt_allele,
                    "wt_population": wt_population,
                    "mut_poulation": mut_poulation, 
                    "wt_freq": wt_population/total_population, # freq should be computed here, since we are decomposing multiple SNVs into separate independent SNVs.
                    "mt_freq": mut_poulation/total_population
                }
                chromosomal_variants.append(v)
        except:
            print(i, tuple)
            raise

        # if i==15:break
    
    return pd.DataFrame(chromosomal_variants)




if __name__ == "__main__":
    print("Log: loading population freq dataset ...")
    inp_filepath = home_dir+"data/ALFA_population_freq/dbsnp_with_prots_and_population_freq_mapping.vcf.gz"
    # for full scale analysis
    dbsnps_df = pd.read_csv(inp_filepath, sep="\t", compression='gzip')

    # for small scale analysis
    small = ""#"_small"
    # dbsnps_iterator = pd.read_csv(inp_filepath, sep="\t", compression='gzip', chunksize=1000)
    # dbsnps_df = dbsnps_iterator.__next__()

    print(dbsnps_df.columns)
    print(f"Log: loaded population freq dataset: {dbsnps_df.shape}")



    colums = ["snp_id", "CHROM", "POS", "REF", "ALT", "SAMN10492705"]
    dbsnps_df = dbsnps_df[colums]
    # print(dbsnps_df.value_counts())


    chromosomal_variants_df = get_chromosomal_variants_df(dbsnps_df)
    chromosomal_variants_df["chrom"] = chromosomal_variants_df["chrom_acc_version"].apply(lambda x: int(x[x.index("_")+1:x.index(".")])) # taking only chromosom number for dbNSFP inputs


    print("Log: saving chromosomal variants w/o population freq in the 'models/aa_common/datasets_population_freq' directory ...")
    chromosomal_variants_df.to_csv(f"models/aa_common/datasets_population_freq/chromosomal_SNVs_with_population_freq{small}.txt", index=False, sep="\t", header=True)
    chromosomal_variants_df[["chrom", "pos", "ref", "alt"]].to_csv(f"models/aa_common/datasets_population_freq/chromosomal_SNVs_without_population_freq{small}.txt", index=False, sep=" ", header=False)
