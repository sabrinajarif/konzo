
BASE = "/scratch/midway3/sjarif/07_MT_exports_HMP"
CONTIGS_DB = "/project/blekhman/shared/konzo/HMP/raw_fastq/MGX/05_CONTIGS_DB"
PROFILE_DB = "/project/blekhman/sjarif/konzo/06_MT_profile_HMP/profiled"
DIAMOND_DB = "/scratch/midway3/sjarif/db/diamond/"


# FIRST: export the gene calls for each contigs db, then filter and extract AA sequences


# export gene calls (works!)
#!/bin/bash -e
#SBATCH -t 36:00:00 -N 1
#SBATCH --mem=32G 
#SBATCH --ntasks=1
#SBATCH -p blekhman
#SBATCH --account pi-blekhman
#SBATCH --mail-type=END
#SBATCH --mail-user=sjarif@uchicago.edu

cd /scratch/midway3/sjarif/07_MT_exports_HMP

anvi-export-gene-calls \
   -c masiManimba.scg.db \
   --gene-caller prodigal \
   -o masiManimba_gene_calls.tsv

# remove gene calls of short sequences (less than 50 AA) and partial sequences
# (works!)
awk -F'\t' '
NR==1 {next}
$6==0 && length($10)>=50 {print $1}
' masiManimba_gene_calls.tsv > masiManimba_gene_ids_min50aa_complete.txt


# get AA sequences for filtered gene calls
# (works!)
#!/bin/bash -e
#SBATCH -t 36:00:00 -N 1
#SBATCH --mem=32G 
#SBATCH --ntasks=1
#SBATCH -p blekhman
#SBATCH --account pi-blekhman
#SBATCH --mail-type=END
#SBATCH --mail-user=sjarif@uchicago.edu

cd /scratch/midway3/sjarif/07_MT_exports_HMP

anvi-get-sequences-for-gene-calls \
    -c masiManimba.scg.db \
    --get-aa-sequences \
    --gene-caller-ids masiManimba_gene_ids_min50aa_complete.txt \
    --output-file masiManimba_genes.faa


# SECOND:  query these AA sequences using diamond. I already have everything realted to dimaond downloaded


# perform the query on the filtered mm contigs
# (works!!!!!!!)

#!/bin/bash -e
#SBATCH -t 36:00:00
#SBATCH -N 1
#SBATCH --mem=75G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH -p caslake
#SBATCH --account pi-blekhman
#SBATCH --mail-type=END
#SBATCH --mail-user=sjarif@uchicago.edu


cd /scratch/midway3/sjarif/db/diamond

diamond blastp \
  --db /scratch/midway3/sjarif/db/diamond/refseq_ref.dmnd \
  --query /scratch/midway3/sjarif/07_MT_exports_HMP/masiManimba_genes.faa \
  --out /scratch/midway3/sjarif/07_MT_exports_HMP/masiManimba_diamond_matches.tsv \
  --outfmt 6 qseqid sseqid pident length evalue bitscore \
  --threads 32 


# THIRD: filter resultant table for best protein accession hits, then merge in taxids


# now filter for the best hit (lowest e-value first, highest bitscore first)
# this does not filter by 90% identity for example, but can do that later in R
# if we get some interesting gene expression patterns we can return to the large tsv
sort -k1,1 -k5,5g -k6,6nr masiManimba_diamond_matches.tsv \
  | awk '!seen[$1]++' \
  > masiManimba_diamond_matches_besthit.tsv

# finally, add taxonomy
# get list of accessions
cut -f2 masiManimba_diamond_matches_besthit.tsv | sort -u > masiManimba_besthit_accessions.txt


zgrep -F -f masiManimba_besthit_accessions.txt ../db/diamond/prot.accession2taxid.gz \
  | awk '{print $1 "\t" $3}' \
  > masiManimba_accession_to_taxid.tsv





#!/bin/bash -e
#SBATCH -t 36:00:00 -N 1
#SBATCH --mem=32G 
#SBATCH --ntasks=1
#SBATCH -p caslake
#SBATCH --account pi-blekhman
#SBATCH --mail-type=END
#SBATCH --mail-user=sjarif@uchicago.edu

set -euo pipefail

BASE="/scratch/midway3/sjarif/07_MT_exports_HMP"
ACC2TAX_DB="/scratch/midway3/sjarif/db/diamond/prot.accession2taxid.gz"

for sample_dir in "$BASE"/*/; do
    sample=$(basename "$sample_dir")

    besthit="${sample_dir}/${sample}_diamond_matches_besthit.tsv"
    accessions="${sample_dir}/${sample}_besthit_accessions.txt"
    out="${sample_dir}/${sample}_accession_to_taxid.tsv"

    # Skip if input is missing
    if [[ ! -f "$besthit" ]]; then
        echo "⚠️  Skipping $sample (missing besthit file)"
        continue
    fi

    # Skip if output already exists
    if [[ -f "$out" ]]; then
        echo "✅ Skipping $sample (already done)"
        continue
    fi

    echo "🔄 Processing $sample"

    cut -f2 "$besthit" | sort -u > "$accessions"

    zgrep -F -f "$accessions" "$ACC2TAX_DB" \
        | awk '{print $1 "\t" $3}' > "$out"
done

echo "🎉 All samples processed"
