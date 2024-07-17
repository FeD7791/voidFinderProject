#!/bin/bash

#SBATCH -o ../Example/popcorn_out.%j
#SBATCH -e ../Example/popcorn_err.%j

#SBATCH -D ../Source/
#SBATCH -J Popcorn

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=72

#SBATCH --time=00:15:00

module purge
module load gcc openmpi hdf5-serial

. /etc/profile
export OMP_NUM_THREADS=72

srun ./svf config=../Example/vars.conf
srun ./popcorn config=../Example/vars.conf
srun ./compute_intersecs config=../Example/vars.conf
srun ./clean_duplicates config=../Example/vars.conf
