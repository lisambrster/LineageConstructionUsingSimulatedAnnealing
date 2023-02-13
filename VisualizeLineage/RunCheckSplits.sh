#!/bin/bash -e

module load matlab/R2021a

#SBATCH --nodes=1
#SBATCH --constraint=skylake
#SBATCH --partition=ccb
matlab -nosplash -nodesktop -r "CheckSplits('CurrentConfig'); exit"



