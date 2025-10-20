#!/bin/bash -e
#SBATCH -t 36:00:00 -N 1
#SBATCH --mem=1G
#SBATCH --ntasks=1
#SBATCH -p blekhman
#SBATCH --account pi-blekhman

# load java
module load java/15.0.2

# set wd
cd /project/blekhman/sjarif/konzo

# activate the virtual environment
#source /project/blekhman/sjarif/konzo_qc/bin/activate
conda init bash
conda activate /project/blekhman/sjarif/conda/envs/konzo_QC

# run snakemake
snakemake --keep-going --jobs 25 --slurm --default-resources slurm_account=pi-blekhman slurm_partition=caslake

# to submit this job:
# sbatch --job-name=next -o nextone.log run_snakemake.slurm
# this will create a log called "nextone.log" in the current directory

