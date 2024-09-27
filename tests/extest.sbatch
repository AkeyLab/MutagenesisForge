#!/bin/bash
#SBATCH --job-name=vepsimtests                           # create a short name for your job
#SBATCH --nodes=1                                # node count
#SBATCH --ntasks=1                               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1                        # cpu-cores per task (equal to $SLURM_CPUS_PER_TASK)
#SBATCH --mem=10G                                # total cpu requested memory
#SBATCH --time=01:00:00                          # total run time limit (HH:MM:SS)
##SBATCH --array=1-10                            # array'd jobs (each job will have $SLURM_ARRAY_TASK_ID)
##SBATCH --output=slurm_outs/slurm-%A_%a.out     # where to store and what to name slurm outs

MutagenesisForge exhaust --fasta /projects/AKEY/akey_vol2/References/Genomes/hs37d5/hs37d5.fa

MutagenesisForge exhaust --fasta /projects/AKEY/akey_vol2/References/Genomes/hs37d5/hs37d5.fa --by-read