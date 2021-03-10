#!/bin/bash
#SBATCH -N 4
#SBATCH -C haswell
#SBATCH -q regular
#SBATCH -J 13f_new_Mixing
#SBATCH --mail-user=fernando_tta@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH -t 12:00:00

module load cray-hdf5/1.10.2.0
#OpenMP settings:
export OMP_NUM_THREADS=64
export OMP_PLACES=threads
export OMP_PROC_BIND=spread


#run the application:
srun -n 4 -c 64 --cpu_bind=cores -C haswell do13f.sh
