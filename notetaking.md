# download multiqc directory
scp -r sjarif@midway3.rcc.uchicago.edu:/project/blekhman/sjarif/konzo/01_MT_initialQC/raw_reads/multiQC /Users/sabri/Documents/GitHub/konzo/preprocessing/01_MT_initialQC


scp -r sjarif@midway3.rcc.uchicago.edu:/project/blekhman/sjarif/konzo/07_MT_exports_HMP/scg_taxonomy_MGX /Users/sabri/Documents/GitHub/konzo/analysis/data



scp -r sjarif@midway3.rcc.uchicago.edu:/project/blekhman/sjarif/konzo/07_MT_exports/gene_coverages/merged_gene_coverages.tsv /Users/sabri/Documents/GitHub/konzo/analysis/data/DRC_MT_diamond_output
scp -r sjarif@midway3.rcc.uchicago.edu:/project/blekhman/sjarif/konzo/07_MT_exports/gene_detection/merged_gene_detection.tsv /Users/sabri/Documents/GitHub/konzo/analysis/data/DRC_MT_diamond_output
scp -r sjarif@midway3.rcc.uchicago.edu:/scratch/midway3/sjarif/07_MT_exports_HMP/masiManimba_diamond_matches_besthit.tsv /Users/sabri/Documents/GitHub/konzo/analysis/data/DRC_MT_diamond_output
scp -r sjarif@midway3.rcc.uchicago.edu:/scratch/midway3/sjarif/07_MT_exports_HMP/masiManimba_accession_to_taxid.tsv /Users/sabri/Documents/GitHub/konzo/analysis/data/DRC_MT_diamond_output

scp -r sjarif@midway3.rcc.uchicago.edu:/scratch/midway3/sjarif/07_MT_exports_HMP/exports /Users/sabri/Documents/GitHub/konzo/analysis/data/HMP_MT_diamond_output





# test snakemake code when not returning logs
snakemake -j 1 --rerun-incomplete --printshellcmds --reason
snakemake -j 1 --rerun-incomplete --printshellcmds --reason --use-conda
snakemake -j 1 --rerun-incomplete --printshellcmds --reason --keep-going

# creating the anvio venv because it was incompatible with my konzo env
cd /project/blekhman/sjarif/konzo
conda create -y -p envs/anvio-8 python=3.10 #envs dir already existed
conda activate anvio-8
conda install -y -c conda-forge -c bioconda python=3.10 \
        sqlite=3.46 prodigal idba mcl muscle=3.8.1551 famsa hmmer diamond \
        blast megahit spades bowtie2 bwa graphviz "samtools>=1.9" \
        trimal iqtree trnascan-se fasttree vmatch r-base r-tidyverse \
        r-optparse r-stringi r-magrittr bioconductor-qvalue meme ghostscript \
        nodejs=20.12.2
conda install -y -c bioconda fastani
curl -L https://github.com/merenlab/anvio/releases/download/v8/anvio-8.tar.gz \
        --output anvio-8.tar.gz
pip install anvio-8.tar.gz
anvi-self-test --version
# had to install a different version of pulp because of errors
pip install pulp==2.7.0



FILENAME_PREFIX=13R_S13

anvi-export-gene-coverage-and-detection \
  -p ${SAMPLE_DIR}/PROFILE.db \
  -c ${CONTIGS_DB} \
  -O ${OUTDIR}/${FILENAME_PREFIX}


  scp -r sjarif@midway3.rcc.uchicago.edu:/project/blekhman/jjcolgan/mgDenovoAssembly/populationAssembly/12_METAGENOME_TAXONOMY/masiManimba /Users/sabri/Documents/GitHub/konzo/preprocessing/06_MT_profile_merge_cluster


echo "[$(date)] Exporting SCG taxonomy table..."

anvi-export-scg-taxonomy \
  -c masiManimba.scg.db \
  --gene-caller prodigal \
  -o scg_taxonomy.txt

echo "[$(date)] SCG taxonomy export finished."
echo "Sanity check: first 5 lines of scg_taxonomy.txt"
head -n 5 scg_taxonomy.txt
echo

echo "[$(date)] Pipeline completed successfully"




### THIS WORKED BUT WITH SOME WARNING FLAGS
anvi-estimate-scg-taxonomy \
  -c masiManimba.scg.db \
  -p ../06_MT_profile_merge_cluster/profiled/mtx/10R_S10/PROFILE.db \
  --compute-scg-coverages \
  --metagenome-mode \
  -o S10R_S10_scg_taxonomy



## I RAN THIS CODE TO LOOP THROUGH ALL THE DRC SAMPLES
CONTIGS_DB="/project/blekhman/sjarif/konzo/07_MT_exports/masiManimba.scg.db"
PROFILE_ROOT="/project/blekhman/sjarif/konzo/06_MT_profile_merge_cluster/profiled/mtx"
OUTDIR="/project/blekhman/sjarif/konzo/07_MT_exports/scg_taxonomy"

mkdir -p "${OUTDIR}"

for SAMPLE_DIR in "${PROFILE_ROOT}"/*; do
    SAMPLE=$(basename "${SAMPLE_DIR}")
    PROFILE_DB="${SAMPLE_DIR}/PROFILE.db"

    # Skip if PROFILE.db does not exist
    if [[ ! -f "${PROFILE_DB}" ]]; then
        echo "Skipping ${SAMPLE}: no PROFILE.db"
        continue
    fi

    echo "Processing sample: ${SAMPLE}"

    anvi-estimate-scg-taxonomy \
        -c "${CONTIGS_DB}" \
        -p "${PROFILE_DB}" \
        --compute-scg-coverages \
        --metagenome-mode \
        -o "${OUTDIR}/${SAMPLE}_scg_taxonomy"

done



## CODE FOR HMP SAMPLES
PROFILE_ROOT="/project/blekhman/sjarif/konzo/06_MT_profile_HMP/profiled"
CONTIGS_ROOT="/project/blekhman/shared/konzo/HMP/raw_fastq/MGX/05_CONTIGS_DB"
OUTDIR="/project/blekhman/sjarif/konzo/07_MT_exports_HMP/scg_taxonomy"

mkdir -p "${OUTDIR}"

for SAMPLE_DIR in "${PROFILE_ROOT}"/*; do
    # Skip non-directories
    if [[ ! -d "${SAMPLE_DIR}" ]]; then
        continue
    fi

    SAMPLE=$(basename "${SAMPLE_DIR}" | xargs)  # trim spaces
    PROFILE_DB="${SAMPLE_DIR}/PROFILE.db"

    # Skip if PROFILE.db does not exist
    if [[ ! -f "${PROFILE_DB}" ]]; then
        echo "Skipping ${SAMPLE}: no PROFILE.db"
        continue
    fi

    # Auto-detect contigs.db in the contigs folder
    CONTIGS_FOLDER="${CONTIGS_ROOT}/${SAMPLE}"
    CONTIGS_DB=$(find "${CONTIGS_FOLDER}" -maxdepth 1 -type f -name "*.db" | head -n 1 || true)

    if [[ -z "${CONTIGS_DB}" ]]; then
        echo "Skipping ${SAMPLE}: no contigs.db found in ${CONTIGS_FOLDER}"
        continue
    fi

    echo "Processing sample: ${SAMPLE}"
    echo "PROFILE_DB = ${PROFILE_DB}"
    echo "CONTIGS_DB  = ${CONTIGS_DB}"

    anvi-estimate-scg-taxonomy \
        -c "${CONTIGS_DB}" \
        -p "${PROFILE_DB}" \
        --compute-scg-coverages \
        --metagenome-mode \
        -o "${OUTDIR}/${SAMPLE}_scg_taxonomy"

done



## CODE FOR HMP MGX SAMPLES
PROFILE_ROOT="/project/blekhman/shared/konzo/HMP/raw_fastq/MGX/06_MG_PROFILES"
CONTIGS_ROOT="/project/blekhman/shared/konzo/HMP/raw_fastq/MGX/05_CONTIGS_DB"
OUTDIR="/project/blekhman/sjarif/konzo/07_MT_exports_HMP/scg_taxonomy_MGX"

mkdir -p "${OUTDIR}"

for SAMPLE_DIR in "${PROFILE_ROOT}"/*; do
    # Skip non-directories
    if [[ ! -d "${SAMPLE_DIR}" ]]; then
        continue
    fi

    SAMPLE=$(basename "${SAMPLE_DIR}" | xargs)  # trim spaces
    PROFILE_DB="${SAMPLE_DIR}/PROFILE.db"

    # Skip if PROFILE.db does not exist
    if [[ ! -f "${PROFILE_DB}" ]]; then
        echo "Skipping ${SAMPLE}: no PROFILE.db"
        continue
    fi

    # Auto-detect contigs.db in the contigs folder
    CONTIGS_FOLDER="${CONTIGS_ROOT}/${SAMPLE}"
    CONTIGS_DB=$(find "${CONTIGS_FOLDER}" -maxdepth 1 -type f -name "*.db" | head -n 1 || true)

    if [[ -z "${CONTIGS_DB}" ]]; then
        echo "Skipping ${SAMPLE}: no contigs.db found in ${CONTIGS_FOLDER}"
        continue
    fi

    echo "Processing sample: ${SAMPLE}"
    echo "PROFILE_DB = ${PROFILE_DB}"
    echo "CONTIGS_DB  = ${CONTIGS_DB}"

    anvi-estimate-scg-taxonomy \
        -c "${CONTIGS_DB}" \
        -p "${PROFILE_DB}" \
        --compute-scg-coverages \
        --metagenome-mode \
        -o "${OUTDIR}/${SAMPLE}_scg_taxonomy"

done

## DRC MGX
CONTIGS_DB="/project/blekhman/sjarif/konzo/07_MT_exports/masiManimba.scg.db"
PROFILE_ROOT="/project/blekhman/sjarif/konzo/06_MT_profile_merge_cluster/profiled/mgx"
OUTDIR="/project/blekhman/sjarif/konzo/07_MT_exports/scg_taxonomy_MGX"

mkdir -p "${OUTDIR}"

for SAMPLE_DIR in "${PROFILE_ROOT}"/*; do
    SAMPLE=$(basename "${SAMPLE_DIR}")
    PROFILE_DB="${SAMPLE_DIR}/PROFILE.db"

    # Skip if PROFILE.db does not exist
    if [[ ! -f "${PROFILE_DB}" ]]; then
        echo "Skipping ${SAMPLE}: no PROFILE.db"
        continue
    fi

    echo "Processing sample: ${SAMPLE}"

    anvi-estimate-scg-taxonomy \
        -c "${CONTIGS_DB}" \
        -p "${PROFILE_DB}" \
        --compute-scg-coverages \
        --metagenome-mode \
        -o "${OUTDIR}/${SAMPLE}_scg_taxonomy"

done

## test for one HMP sample

cd /scratch/midway3/sjarif/CSM79HGX_test

anvi-export-gene-coverage-and-detection \
  -c /project/blekhman/shared/konzo/HMP/raw_fastq/MGX/05_CONTIGS_DB/CSM79HGX/contigs.db \
  -p /project/blekhman/sjarif/konzo/06_MT_profile_HMP/profiled/CSM79HGX/PROFILE.db \
  -O CSM79HGX

anvi-get-sequences-for-gene-calls \
  -c /project/blekhman/shared/konzo/HMP/raw_fastq/MGX/05_CONTIGS_DB/CSM79HGX/contigs.db \
  -o gene_calls.fa

kaiju-makedb -s refseq_ref


anvi-export-contigs \
  -c /project/blekhman/shared/konzo/HMP/raw_fastq/MGX/05_CONTIGS_DB/CSM79HGX/contigs.db \
  -o CSM79HGX_taxonomy





# metaphlan/ humann

conda activate /scratch/midway3/sjarif/metaphlan_venv

humann_databases \
  --download chocophlan full /scratch/midway3/sjarif/db/humann \
  --update-config yes 

humann_databases \
  --download uniref uniref90_diamond /scratch/midway3/sjarif/db/humann \
  --update-config yes


# concat reads

#!/bin/bash
#SBATCH --job-name=mg_concat_fastq
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH -p caslake
#SBATCH --account pi-blekhman


# INPUT directory with paired-end FASTQs
IN_DIR="/project/blekhman/shared/konzo/HMP/raw_fastq/MGX/01_QC"

# OUTPUT directory for concatenated FASTQs
OUT_DIR="/scratch/midway3/sjarif/MG_concat_fastq_HMP"


# Loop over R1 files
for R1 in ${IN_DIR}/*_R1*_dehosted.fastq.gz; do

    # Infer R2 filename
    R2="${R1/_R1/_R2}"

    # Extract sample name
    SAMPLE=$(basename "${R1}" | sed 's/_R1.*//')

    # Output file
    OUT_FASTQ="${OUT_DIR}/${SAMPLE}.fastq.gz"

    # Check that R2 exists
    if [[ ! -f "${R2}" ]]; then
        echo "WARNING: R2 file not found for ${SAMPLE}, skipping"
        continue
    fi

    echo "Concatenating ${SAMPLE}"

    # Concatenate while preserving gzip
    zcat "${R1}" "${R2}" | gzip > "${OUT_FASTQ}"

done













# snakefile for metaphlan on all hmp samples
# THIS WORKS

import os
MGX_DIR = "/scratch/midway3/sjarif/MG_concat_fastq_HMP"
OUT_DIR = "/scratch/midway3/sjarif/08_MG_metaphlan_HMP"
METAPHLAN_DB = "/scratch/midway3/sjarif/db/metaphlan"

SAMPLES = [os.path.basename(f).replace(".fastq.gz", "") 
           for f in os.listdir(MGX_DIR) if f.endswith(".fastq.gz")]


rule all:
    input:
        expand(os.path.join(OUT_DIR, "{sample}_out.tsv"), sample=SAMPLES)


rule run_metaphlan:
    input:
        fastq=os.path.join(MGX_DIR, "{sample}.fastq.gz")
    output:
        profile=os.path.join(OUT_DIR, "{sample}_out.tsv")
    params:
        db=METAPHLAN_DB
    threads: 8
    resources:
        mem_mb=64000
    shell:
        """
        metaphlan {input.fastq} \
            --input_type fastq \
            --bowtie2db {params.db} \
            --index mpa_v30_CHOCOPhlAn_201901 \
            --nproc {threads} \
            -o {output.profile}
        """



## run metaphlan on one sample

metaphlan /scratch/midway3/sjarif/MG_concat_fastq_HMP/CSM79HGX.fastq.gz \
    --input_type fastq \
    --bowtie2db /scratch/midway3/sjarif/db/metaphlan \
    --index mpa_v30_CHOCOPhlAn_201901 \
    --nproc 8 \
    -o /scratch/midway3/sjarif/08_MG_metaphlan_HMP/CSM79HGX_out.tsv


## run humann on one mtx hmp sample
## I THINK THIS WORKS
humann \
  --input /scratch/midway3/sjarif/MT_concat_fastq_HMP/CSM79HGX.fastq.gz \
  --input-format fastq.gz \
  --taxonomic-profile /scratch/midway3/sjarif/08_MG_metaphlan_HMP/CSM79HGX_out.tsv \
  --nucleotide-database /scratch/midway3/sjarif/db/humann/chocophlan \
  --output /scratch/midway3/sjarif/09_MT_humann_HMP/CSM79HGX \
  --threads 8



# snakefile for humann on all hmp MTX samples

import os
MTX_DIR = "/scratch/midway3/sjarif/MT_concat_fastq_HMP"
OUT_DIR = "/scratch/midway3/sjarif/09_MT_humann_HMP"
TAX_PROFILE = "/scratch/midway3/sjarif/08_MG_metaphlan_HMP"
NT_DB = "/scratch/midway3/sjarif/db/humann/chocophlan"

SAMPLES = [os.path.basename(f).replace("_out.tsv", "") 
           for f in os.listdir(TAX_PROFILE) if f.endswith("_out.tsv")]


rule all:
    input:
        expand(os.path.join(OUT_DIR, "{sample}", "{sample}_genefamilies.tsv"), sample=SAMPLES)


rule run_humann:
    input:
        fastq=os.path.join(MTX_DIR, "{sample}.fastq.gz"),
        taxprofile=os.path.join(TAX_PROFILE, "{sample}_out.tsv")
    output:
        genefam=os.path.join(OUT_DIR, "{sample}", "{sample}_genefamilies.tsv")
    params:
        outdir=os.path.join(OUT_DIR, "{sample}"),
        db=NT_DB
    threads: 8
    resources:
        mem_mb=64000
    shell:
        """
        humann \
            --input {input.fastq} \
            --input-format fastq.gz \
            --taxonomic-profile {input.taxprofile} \
            --nucleotide-database {params.db} \
            --output {params.outdir} \
            --threads {threads}
        """


## humann code to join all tables together

humann_join_tables \
  --input /scratch/midway3/sjarif/09_MT_humann_HMP \
  --file_name genefamilies \
  --search-subdirectories \
  --output humann_all_genefamilies_HMP.tsv

humann_join_tables \
  --input /scratch/midway3/sjarif/09_MT_humann_HMP \
  --file_name pathabundance \
  --search-subdirectories \
  --output humann_all_pathabundance_HMP.tsv

humann_join_tables \
  --input /scratch/midway3/sjarif/09_MT_humann_HMP \
  --file_name pathcoverage \
  --search-subdirectories \
  --output humann_all_pathcoverage_HMP.tsv





snakemake -j 1 --rerun-incomplete --printshellcmds --reason








# check hmp contigs
ssh -L 8080:localhost:8080 sjarif@midway3.rcc.uchicago.edu

cd /project/blekhman/shared/konzo/HMP/raw_fastq/MGX/05_CONTIGS_DB
anvi-display-contigs-stats */contigs.db

