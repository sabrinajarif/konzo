#!/bin/bash
#SBATCH --job-name=anvio-db
#SBATCH --output=anvio-db-%j.log
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -p blekhman
#SBATCH --account pi-blekhman
#SBATCH -t 36:00:00 -N 1
#SBATCH --ntasks=1

cd /project/blekhman/sjarif/konzo
conda activate anvio-8

# Define input/output
FASTA="/project/blekhman/shared/konzo/contigs/masiManimba.contigs.simplified.fa"
DB_OUT="/project/blekhman/shared/konzo/contigs/masiManimba.db"
NAME="masiManimba contigs"

# Run anvi-gen-contigs-database
anvi-gen-contigs-database -f "$FASTA" -o "$DB_OUT" -n "$NAME" --num-threads $SLURM_CPUS_PER_TASK
