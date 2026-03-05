#!/bin/bash -e
#SBATCH -t 36:00:00 -N 1
#SBATCH --mem=1G
#SBATCH --ntasks=1
#SBATCH -p caslake
#SBATCH --account pi-blekhman
#SBATCH --mail-type=END
#SBATCH --mail-user=sjarif@uchicago.edu

# run snakemake
snakemake \
  --keep-going \
  --rerun-incomplete \
  --jobs 20 \
  --slurm \
  --default-resources \
    slurm_account=pi-blekhman \
    slurm_partition=caslake \
    mem_mb=2000


# to submit this job:
# navigate to cd /project/blekhman/sjarif/konzo
# conda activate anvio-dev
# cd 06_MT_profile_HMP
# sbatch --job-name=06_MT_profile_HMP -o 006_MT_profile_HMP.log run_snakemake.slurm
# this will create a log called "06_MT_profile_HMP.log" in the current directory