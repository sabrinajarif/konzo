#!/bin/bash -e
#SBATCH -t 36:00:00 -N 1
#SBATCH --mem=48G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH -p caslake
#SBATCH --account pi-blekhman
#SBATCH --mail-type=END
#SBATCH --mail-user=sjarif@uchicago.edu

source ~/.bashrc
cd /project/blekhman/sjarif/konzo
conda activate anvio-dev
cd 07_MT_exports

echo "Running HMMs..."

anvi-run-hmms \
  -c masiManimba.scg.db 

echo "HMMs finished."
