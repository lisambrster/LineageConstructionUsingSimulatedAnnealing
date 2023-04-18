#!/bin/bash
#SBATCH --nodes=1
#sbatch -p ccb run_MakeRegisteredImages.sh
#SBATCH --ntasks-per-node=1
#SBATCH -t 7-00:00            # wall time (D-HH:MM)
#SBATCH --output=logs/log_MakeRegisteredImages.out

export HDF5_USE_FILE_LOCKING='FALSE'

python3 MakeRegisteredImages.py -c config.yaml
