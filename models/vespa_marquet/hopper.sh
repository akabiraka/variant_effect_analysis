#!/bin/bash
#SBATCH --job-name=VESPA
#SBATCH --output=/projects/ashehu/akabir4/projects/variant_effect_analysis/outputs/logs/popu_freq-%j.out
#SBATCH --error=/projects/ashehu/akabir4/projects/variant_effect_analysis/outputs/logs/popu_freq-%j.err
#SBATCH --mail-user=<akabir4@gmu.edu>
#SBATCH --mail-type=BEGIN,END,FAIL

## cpu 
##SBATCH --partition=contrib                  # submit   to the normal(default) partition
##SBATCH --mem=16GB                # Request nGB RAM per core
##SBATCH --time=0-08:00   # Total time needed for job: Days-Hours:Minutes

#SBATCH --partition=gpuq 
#SBATCH --qos=gpu 
#SBATCH --gres=gpu:A100.40gb:1 
#SBATCH --mem=64G

source /projects/ashehu/akabir4/venvs/hopper_vespa_marquet/bin/activate

vespa_logodds models/vespa_marquet/cache/merged_popu_and_patho_sequences.fasta --prott5_weights_cache models/vespa_marquet/cache -o models/vespa_marquet/logodds/popu_and_patho_logodds.h5 --csv_dir models/vespa_marquet/logodds_popu_and_patho/

##vespa_logodds models/aa_common/datasets_pmd_analysis/pmd_sequences.fasta --prott5_weights_cache models/vespa_marquet/cache -o models/vespa_marquet/logodds/pmd_logodds.h5 --csv_dir models/vespa_marquet/logodds_pmd/