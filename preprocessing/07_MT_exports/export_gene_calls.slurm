#!/bin/bash -e
#SBATCH -t 36:00:00 -N 1
#SBATCH --mem=100G
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


echo "Exporting gene calls (prodigal)..."

anvi-export-gene-calls \
  -c masiManimba.scg.db \
  --gene-caller prodigal \
  -o gene_calls.txt

echo "Gene calls export finished."
echo "Sanity check: first 3 lines of gene_calls.txt"
head -n 3 gene_calls.txt
echo

