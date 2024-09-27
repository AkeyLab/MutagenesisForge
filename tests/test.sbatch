#!/bin/bash
#SBATCH --job-name=vepsimtests                           # create a short name for your job
#SBATCH --nodes=1                                # node count
#SBATCH --ntasks=1                               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1                        # cpu-cores per task (equal to $SLURM_CPUS_PER_TASK)
#SBATCH --mem=10G                                # total cpu requested memory
#SBATCH --time=10:00:00                          # total run time limit (HH:MM:SS)
##SBATCH --array=1-10                            # array'd jobs (each job will have $SLURM_ARRAY_TASK_ID)
##SBATCH --output=slurm_outs/slurm-%A_%a.out     # where to store and what to name slurm outs

MutagenesisForge context --vcf /projects/AKEY/akey_vol2/huixinx/Projects/01.eGTEx/NWGC/04.fig3/02.exp_mis_to_syn_ratio/docker_stringent.nwgc.rep2.raw_bb_p_lt_10_8.filtered10.txt --bed /projects/AKEY/akey_vol2/huixinx/Projects/01.eGTEx/NWGC/04.fig3/02.exp_mis_to_syn_ratio/step12.problematic.bed --fasta /projects/AKEY/akey_vol2/References/Genomes/hs37d5/hs37d5.fa --output tester --sims 10 --vep_call False
