#!/bin/bash -e
#SBATCH -t 36:00:00 -N 1
#SBATCH --mem=1G
#SBATCH --ntasks=1
#SBATCH -p caslake
#SBATCH --account pi-blekhman
#SBATCH --mail-type=END
#SBATCH --mail-user=sjarif@uchicago.edu

# run snakemake
snakemake --keep-going --rerun-incomplete --jobs 20 --slurm --default-resources slurm_account=pi-blekhman slurm_partition=caslake

# to submit this job:
# navigate to cd /project/blekhman/sjarif/konzo
# conda activate anvio-dev
# cd /scratch/midway3/sjarif/07_MT_exports_HMP
# sbatch --job-name=07_MT_exports_HMP -o 07_MT_exports_HMP.log run_snakemake.slurm
