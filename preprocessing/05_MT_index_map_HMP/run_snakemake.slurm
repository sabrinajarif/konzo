#!/bin/bash -e
#SBATCH -t 36:00:00 -N 1
#SBATCH --mem=1G
#SBATCH --ntasks=1
#SBATCH -p caslake
#SBATCH --account pi-blekhman
#SBATCH --mail-type=END
#SBATCH --mail-user=sjarif@uchicago.edu

# run snakemake
snakemake --keep-going --rerun-incomplete --jobs 20 --slurm --default-resources slurm_account=pi-blekhman slurm_partition=blekhman

# to submit this job:
# navigate to cd /project/blekhman/sjarif/konzo
# conda activate konzo
# cd 05_MT_map_profile_HMP
# sbatch --job-name=05_MT_index_map_HMP -o 05_MT_index_map_HMP.log run_snakemake.slurm
# this will create a log called "05_MT_index_map_HMP.log" in the current directory

