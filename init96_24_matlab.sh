#!/bin/bash
#SBATCH --partition=mammoth         # the requested queue
#SBATCH --nodelist=f03-24
#SBATCH --job-name=THALAMOCORTICAL  # name the job           
#SBATCH --nodes=1                   # number of nodes to use
#SBATCH --ntasks=96                 # total number of tasks (preocesses)
#SBATCH --mem-per-cpu=12000         # in megabytes, unless unit explicitly stated
#SBATCH --error=%J.err              # redirect stderr to this file
#SBATCH --output=%J.out             # redirect stdout to this file
#SBATCH --mail-user=martynas.dervinis@gmail.com     # email address used for event notification
##SBATCH --mail-type=start                          # email on job start  
#SBATCH --mail-type=end                             # email on job end
#SBATCH --mail-type=fail                            # email on job failure
#SBATCH --time=4800

# module load  neuron/7.4

echo "Usable Environment Variables:"
echo "============================="
echo \$SLURM_JOB_ID=${SLURM_JOB_ID} 
echo \$SLURM_NTASKS=${SLURM_NTASKS}
echo \$SLURM_NTASKS_PER_NODE=${SLURM_NTASKS_PER_NODE}
echo \$SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}
echo \$SLURM_JOB_CPUS_PER_NODE=${SLURM_JOB_CPUS_PER_NODE}
echo \$SLURM_MEM_PER_CPU=${SLURM_MEM_PER_CPU}
echo "module list:"
module list 2>&1

# Some of these environment variables are utilised by the qe executable itself
export MYJOB_DATAPATH=${PWD}

# execute the parallel job (we also time it)
time matlab < AnPSD.m