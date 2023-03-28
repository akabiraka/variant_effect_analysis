#!/bin/bash
#SBATCH --job-name=dbNSFP
#SBATCH --output=/projects/ashehu/akabir4/projects/variant_effect_analysis/outputs/logs/popu_freq-%j.out
#SBATCH --error=/projects/ashehu/akabir4/projects/variant_effect_analysis/outputs/logs/popu_freq-%j.err
#SBATCH --mail-user=<akabir4@gmu.edu>
#SBATCH --mail-type=BEGIN,END,FAIL

## cpu 
#SBATCH --partition=contrib                  # submit   to the normal(default) partition
#SBATCH --mem=16GB                # Request nGB RAM per core
#SBATCH --time=0-08:00   # Total time needed for job: Days-Hours:Minutes


## Load the relevant modules needed for the job
module load openjdk # openjdk/11.0.2-qg
source /projects/ashehu/akabir4/venvs/hopper_variant_effect_analysis_mine/bin/activate

## for population frequency analysis
##python models/aa_common/chromosomal_SNV_conversion_for_population_freq.py

## for full scale run it took ~32 minutes
# java -cp data/dbnsfp/dbNSFP43a/ search_dbNSFP43a -i models/dbnsfp/datasets_population_freq/SNVs_with_popu_freq_balanced_chromosomal.txt -o models/dbnsfp/outputs/popu_freq_preds.txt -w 1-6,77,90,38,44,47,118,83 


## for pathogenicity analysis
for i in {7..9}
do
    # java -cp data/dbnsfp/dbNSFP43a/ search_dbNSFP43a -i models/dbnsfp/datasets_pathogenicity/pathogenic_and_neutral_SNVs/$i.txt -o models/dbnsfp/outputs/pathogenic_and_neutral_SNVs/$i.txt -w 1-6,77,90,38,44,47,118,83
    java -cp data/dbnsfp/dbNSFP43a/ search_dbNSFP43a -i models/dbnsfp/datasets_pathogenicity/likely_pathogenic_and_neutral_SNVs/$i.txt -o models/dbnsfp/outputs/likely_pathogenic_and_neutral_SNVs/$i.txt -w 1-6,77,90,38,44,47,118,83
done