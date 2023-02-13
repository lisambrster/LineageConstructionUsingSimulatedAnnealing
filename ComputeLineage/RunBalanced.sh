#!/bin/bash -e

id=$1
shift

module load slurm  python  disBatch/beta
#sbatch -p ccb -n 6 --ntasks-per-node 1 disBatch --logfile ./disBatchLogs/log.txt taskfile.txt

python AnnealBalanced.py  "$@" 2> logs/${id}.log