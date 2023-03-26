#!/bin/bash
#SBATCH --job-name=dbNSFP
#SBATCH --output=/projects/ashehu/akabir4/projects/variant_effect_analysis/outputs/logs/popu_freq-%j.out
#SBATCH --error=/projects/ashehu/akabir4/projects/variant_effect_analysis/outputs/logs/popu_freq-%j.err
#SBATCH --mail-user=<akabir4@gmu.edu>
#SBATCH --mail-type=BEGIN,END,FAIL

## cpu 
#SBATCH --partition=normal                  # submit   to the normal(default) partition
#SBATCH --mem=32GB                # Request nGB RAM per core
#SBATCH --time=0-02:00   # Total time needed for job: Days-Hours:Minutes


## Load the relevant modules needed for the job
module load openjdk # openjdk/11.0.2-qg
source /projects/ashehu/akabir4/venvs/hopper_variant_effect_analysis_mine/bin/activate

## for population frequency analysis
##python models/aa_common/chromosomal_SNV_conversion_for_population_freq.py

## for full scale run it took ~32 minutes
java -cp data/dbnsfp/dbNSFP43a/ search_dbNSFP43a -i models/dbnsfp/inputs/chromosomal_SNVs.txt -o models/dbnsfp/outputs/chromosomal_SNVs_preds.txt -w 1-6,77-79,90-91
##python models/dbnsfp/analysis_population_freq.py


## for pathogenicity analysis
## java -cp data/dbnsfp/dbNSFP43a/ search_dbNSFP43a -i models/aa_common/inputs/clinvar_HumanPathogenicMissenseVariants01012022To14022023.txt -o models/dbnsfp/outputs/clinvar_HumanPathogenicMissenseVariants01012022To14022023.out -w 1-6,77-79,85-86,90-91,94-96,124-125    
## java -cp data/dbnsfp/dbNSFP43a/ search_dbNSFP43a -i models/aa_common/inputs/clinvar_HumanLikelyPathogenicMissenseVariants01012022To14022023.txt -o models/dbnsfp/outputs/clinvar_HumanLikelyPathogenicMissenseVariants01012022To14022023.out -w 1-6,77-79,85-86,91-92,94-96,124-125