#!/bin/bash -e
#SBATCH -t 36:00:00 -N 1
#SBATCH --mem=1G
#SBATCH --ntasks=1
#SBATCH -p blekhman
#SBATCH --account pi-blekhman
#SBATCH --mail-type=END
#SBATCH --mail-user=sjarif@uchicago.edu

# run snakemake
snakemake --keep-going --jobs 20 --slurm --default-resources slurm_account=pi-blekhman slurm_partition=blekhman

# to submit this job:
# navigate to cd /project/blekhman/sjarif/konzo
# conda activate konzo
# cd 02_MT_trim_dehost_rRNAremove_HMP
# sbatch --job-name=02_MT_trim_dehost_rRNAremove_HMP -o 02_MT_trim_dehost_rRNAremove_HMP.log run_snakemake.slurm
# this will create a log called "02_MT_trim_dehost_rRNAremove_HMP.log" in the current directory

