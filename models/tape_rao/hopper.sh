#!/bin/bash
#SBATCH --job-name=protbert
#SBATCH --output=/projects/ashehu/akabir4/projects/variant_effect_analysis/outputs/logs/popu_freq-%j.out
#SBATCH --error=/projects/ashehu/akabir4/projects/variant_effect_analysis/outputs/logs/popu_freq-%j.err
#SBATCH --mail-user=<akabir4@gmu.edu>
#SBATCH --mail-type=BEGIN,END,FAIL

## cpu 
#SBATCH --partition=normal                  # submit   to the normal(default) partition
#SBATCH --mem-per-cpu=16G
#SBATCH --ntasks=11               # Request nGB RAM per core


## Load the relevant modules needed for the job
module load gnu10
module load openmpi
source /projects/ashehu/akabir4/venvs/hopper_tape_rao/bin/activate
mpirun -np 9 python -m mpi4py.futures models/tape_rao/popu_freq_pred.py