#!/bin/bash
#
#PBS -N M191128_B_MD
#PBS -q parallel
#PBS -l walltime=48:00:00
#PBS -l vmem=60gb
#PBS -l nodes=1:ppn=28
#PBS -m bea
#PBS -M md406@le.ac.uk

# Set OMP_NUM_THREADS for OpenMP jobs
export OMP_NUM_THREADS=$PBS_NUM_PPN

# Execute the job code
module load matlab/2018a
cd /scratch/neuronal/md406/runNSG_M191128_B_MD
matlab -nodisplay -nojvm -r AnPSD