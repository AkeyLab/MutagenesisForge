#!/bin/bash
#SBATCH --job-name=exhausttest                           # create a short name for your job
#SBATCH --nodes=1                                # node count
#SBATCH --ntasks=1                               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1                        # cpu-cores per task (equal to $SLURM_CPUS_PER_TASK)
#SBATCH --mem=10G                                # total cpu requested memory
#SBATCH --time=10:00:00                          # total run time limit (HH:MM:SS)
#SBATCH --output=slurm-%x-%j.out                 # for non-array job, name slurm-out with job-name (%x) and job-num (%j)
##SBATCH --array=1-10                            # array'd jobs (each job will have $SLURM_ARRAY_TASK_ID)
##SBATCH --output=slurm_outs/slurm-%x-%A_%a.out  # for array job, name slurm-out with job-name (%x), job-id (%A), and array-num (%a)

#set bash strict mode (http://redsymbol.net/articles/unofficial-bash-strict-mode/)
set -euo pipefail

MutagenesisForge exhaust --fasta chr21.fa
