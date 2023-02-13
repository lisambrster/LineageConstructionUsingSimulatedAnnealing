#!/bin/bash -e

module load matlab/R2021a


cd /mnt/home/lbrown/LineageConstruction/CPD2/core
matlab -nosplash -nodesktop -r "cpd_make( ); exit"
cd /mnt/home/lbrown/LineageConstruction/VisualizeLineage/


#SBATCH --nodes=1
#SBATCH --constraint=skylake
#SBATCH --partition=ccb
matlab -nosplash -nodesktop -r "MakeGraphAndCheckMatches('CurrentConfig'); exit"
#matlab -nosplash -nodesktop -r "VisualizeSequence_LB( ); exit"
# to run without script, just load matlab and run last line..


