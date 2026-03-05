#!/bin/bash -e
#SBATCH -t 36:00:00 -N 1
#SBATCH --mem=1G
#SBATCH --ntasks=1
#SBATCH -p blekhman
#SBATCH --account pi-blekhman
#SBATCH --mail-type=END
#SBATCH --mail-user=sjarif@uchicago.edu

# run snakemake
snakemake --keep-going --jobs 25 --slurm --default-resources slurm_account=pi-blekhman slurm_partition=blekhman

# to submit this job:
# navigate to cd /project/blekhman/sjarif/konzo
# conda activate konzo
# cd 01_MGMT_downloadHMP
# sbatch --job-name=01_MGMT_downloadHMP -o 01_MGMT_downloadHMP run_snakemake.slurm
# this will create a log called "01_MGMT_downloadHMP.log" in the current directory