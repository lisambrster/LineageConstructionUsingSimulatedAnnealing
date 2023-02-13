#!/bin/bash -e

module load matlab/R2021a  #


#SBATCH --nodes=1
#SBATCH --constraint=skylake
#SBATCH --partition=ccb
matlab -nosplash -nodesktop -r "MakeFeatures('CurrentConfig'); exit"
# to run without script, just load matlab and run last line..



