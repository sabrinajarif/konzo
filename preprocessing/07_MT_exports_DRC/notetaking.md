# BASE = "/project/blekhman/sjarif/konzo/07_MT_exports_HMP"
CONTIGS_DB = "/project/blekhman/shared/konzo/HMP/raw_fastq/MGX/05_CONTIGS_DB"
PROFILE_BASE = "/project/blekhman/sjarif/konzo/06_MT_profile_HMP/profiled"

# confirm contigs.db exists and contains useful info, for one sample

anvi-db-info /project/blekhman/shared/konzo/HMP/raw_fastq/MGX/05_CONTIGS_DB/MSM79H54/contigs.db
                                                                                                                                                       
DB Info (no touch)
===============================================
Database Path ................................: /project/blekhman/shared/konzo/HMP/raw_fastq/MGX/05_CONTIGS_DB/MSM79H54/contigs.db                     
description ..................................: [Not found, but it's OK]                                                                               
db_type ......................................: contigs (variant: unknown)                                                                             
version ......................................: 24                                                                                                     

                                                                                                                                                       
DB Info (no touch also)
===============================================
project_name .................................: contigs_simplified                                                                                     
contigs_db_hash ..............................: hash01bed99f                                                                                           
split_length .................................: 20000                                                                                                  
kmer_size ....................................: 4                                                                                                      
num_contigs ..................................: 23028                                                                                                  
total_length .................................: 56764998                                                                                               
num_splits ...................................: 23493                                                                                                  
gene_level_taxonomy_source ...................: None                                                                                                   
genes_are_called .............................: 1                                                                                                      
external_gene_calls ..........................: 0                                                                                                      
external_gene_amino_acid_seqs ................: 0                                                                                                      
skip_predict_frame ...........................: 0                                                                                                      
splits_consider_gene_calls ...................: 1                                                                                                      
trna_taxonomy_was_run ........................: 0                                                                                                      
trna_taxonomy_database_version ...............: None                                                                                                   
reaction_network_ko_annotations_hash .........: None                                                                                                   
reaction_network_kegg_database_release .......: None                                                                                                   
reaction_network_modelseed_database_sha ......: None                                                                                                   
reaction_network_consensus_threshold .........: None                                                                                                   
reaction_network_discard_ties ................: None                                                                                                   
creation_date ................................: 1763135616.26542                                                                                       
gene_function_sources ........................: KEGG_BRITE,KEGG_Module,Pfam,KEGG_Class,KOfam                                                           
modules_db_hash ..............................: d20a0dcd2128                                                                                           
scg_taxonomy_was_run .........................: 1                                                                                                      
scg_taxonomy_database_version ................: GTDB: v214.1; Anvi'o: v1                                                                               
                                                                                                                                                       
* Please remember that it is never a good idea to change these values. But in some
  cases it may be absolutely necessary to update something here, and a
  programmer may ask you to run this program and do it. But even then, you
  should be extremely careful.

                                                                                                                                                       
AVAILABLE GENE CALLERS
===============================================
* 'prodigal' (66,673 gene calls)                                                                                                                       
* 'Ribosomal_RNA_23S' (29 gene calls)                                                                                                                  
* 'Ribosomal_RNA_16S' (12 gene calls)                                                                                                                  

                                                                                                                                                       
AVAILABLE FUNCTIONAL ANNOTATION SOURCES
===============================================
* KEGG_BRITE (24,285 annotations)                                                                                                                      
* KEGG_Class (4,790 annotations)                                                                                                                       
* KEGG_Module (4,790 annotations)                                                                                                                      
* KOfam (24,317 annotations)                                                                                                                           
* Pfam (84,773 annotations)                                                                                                                            

                                                                                                                                                       
AVAILABLE HMM SOURCES
===============================================
* 'Archaea_76' (76 models with 624 hits)                                                                                                               
* 'Bacteria_71' (71 models with 1,162 hits)                                                                                                            
* 'Protista_83' (83 models with 45 hits)                                                                                                               
* 'Ribosomal_RNA_12S' (1 model with 0 hits)                                                                                                            
* 'Ribosomal_RNA_16S' (3 models with 12 hits)                                                                                                          
* 'Ribosomal_RNA_18S' (1 model with 0 hits)                                                                                                            
* 'Ribosomal_RNA_23S' (2 models with 29 hits)                                                                                                          
* 'Ribosomal_RNA_28S' (1 model with 0 hits)                                                                                                            
* 'Ribosomal_RNA_5S' (5 models with 0 hits)             




# so i don't have gene-level taxonomy annotations. Let's do this with kaiju

# run this in /project/blekhman/sjarif/konzo 
conda activate anvio-dev

# copy contigs over to scratch (works!)
cp /project/blekhman/shared/konzo/HMP/raw_fastq/MGX/05_CONTIGS_DB/MSM79H54/contigs.db /scratch/midway3/sjarif/07_MT_exports_HMP/MSM79H54_contigs.db

cp /project/blekhman/sjarif/konzo/07_MT_exports/masiManimba.scg.db /scratch/midway3/sjarif/07_MT_exports_HMP/masiManimba.scg.db

# Example for sample MSM79H54. Get gene AA seqs (works!)
anvi-get-sequences-for-gene-calls \
    -c /scratch/midway3/sjarif/07_MT_exports_HMP/MSM79H54_contigs.db \
    --get-aa-sequences \
    -o /scratch/midway3/sjarif/07_MT_exports_HMP/MSM79H54_genes.faa

anvi-get-sequences-for-gene-calls \
    -c /scratch/midway3/sjarif/07_MT_exports_HMP/masiManimba.scg.db \
    --get-aa-sequences \
    -o /scratch/midway3/sjarif/07_MT_exports_HMP/masiManimba_genes.faa

# install kaiju (works!)
conda install kaiju
kaiju --version
kaiju: invalid option -- '-'
Kaiju 1.10.1
Copyright 2015-2023 Peter Menzel, Anders Krogh
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>

Usage:
   kaiju -t nodes.dmp -f kaiju_db.fmi -i reads.fastq [-j reads2.fastq]

Mandatory arguments:
   -t FILENAME   Name of nodes.dmp file
   -f FILENAME   Name of database (.fmi) file
   -i FILENAME   Name of input file containing reads in FASTA or FASTQ format

Optional arguments:
   -j FILENAME   Name of second input file for paired-end reads
   -o FILENAME   Name of output file. If not specified, output will be printed to STDOUT
   -z INT        Number of parallel threads for classification (default: 1)
   -a STRING     Run mode, either "mem"  or "greedy" (default: greedy)
   -e INT        Number of mismatches allowed in Greedy mode (default: 3)
   -m INT        Minimum match length (default: 11)
   -s INT        Minimum match score in Greedy mode (default: 65)
   -E FLOAT      Minimum E-value in Greedy mode (default: 0.01)
   -x            Enable SEG low complexity filter (enabled by default)
   -X            Disable SEG low complexity filter
   -p            Input sequences are protein sequences
   -v            Enable verbose output

# Download kaiju db in a slurm job with lots of memory. Only need to do this once hopefully
# Originally I chose refseq_nr to ensure detection of rare taxa in the gene expression data
# However even on a 220G 4 thread job it was running out of memory
# Switching to progenoems instead. URL was outdated so i manually downlaoded it
# (works!)

mkdir /scratch/midway3/sjarif/07_MT_exports_HMP/kaijudb
cd /scratch/midway3/sjarif/07_MT_exports_HMP/kaijudb

wget https://progenomes.embl.de/data/repGenomes/pg4_proteins_representatives.faa.gz
gunzip pg4_proteins_representatives.faa.gz

ls -lh pg4_proteins_representatives.faa
head -n 2 pg4_proteins_representatives.faa


# rename so its compatible with kaiju
# (works!)
mv kaiju_db_pg4_proteins_representatives.faa kaiju_db_progenomes.faa
mv kaiju_db_progenomes.faa progenomes


# slurm file to build index
# (works!)
#!/bin/bash -e
#SBATCH -t 36:00:00 -N 1
#SBATCH --mem=220G 
#SBATCH --ntasks=1
#SBATCH -p blekhman
#SBATCH --account pi-blekhman
#SBATCH --mail-type=END
#SBATCH --mail-user=sjarif@uchicago.edu

cd /scratch/midway3/sjarif/07_MT_exports_HMP/kaijudb

kaiju-makedb \
  -s progenomes \
  --index-only 




# Check that nodes.dmp and kaiju_db.fmi exist
# (works!)
ls /scratch/midway3/sjarif/07_MT_exports_HMP/kaijudb
ls -lh /scratch/midway3/sjarif/07_MT_exports_HMP/kaijudb/progenomes | grep -E "bwt|sa"



# Run kaiju HMP
#!/bin/bash -e
#SBATCH -t 36:00:00 -N 1
#SBATCH --mem=150G 
#SBATCH --ntasks=1
#SBATCH -p blekhman
#SBATCH --account pi-blekhman
#SBATCH --mail-type=END
#SBATCH --mail-user=sjarif@uchicago.edu

kaiju -t /scratch/midway3/sjarif/07_MT_exports_HMP/kaijudb/nodes.dmp \
      -f /scratch/midway3/sjarif/07_MT_exports_HMP/kaijudb/progenomes/kaiju_db_progenomes.fmi \
      -i /scratch/midway3/sjarif/07_MT_exports_HMP/MSM79H54_genes.faa \
      -o /scratch/midway3/sjarif/07_MT_exports_HMP/MSM79H54_kaiju.out \
      -z 16 \
      -a greedy \
      -e 5 \
      -p \
      -v

# add names
kaiju-addTaxonNames -t kaijudb/nodes.dmp -n kaijudb/names.dmp -i MSM79H54_kaiju.out -o MSM79H54_kaiju-names.out

kaiju-addTaxonNames -t ../db/kaijudb/nodes.dmp -n ../db/kaijudb/names.dmp -i masiManimba_kaiju.out -o masiManimba_kaiju-names.out


#### this returned ALL unclassified genes. No idea why
#### I tried adding taxonomy (below), but that doesnt work for unclassified
#### swtiching to MM contigs db because it was a co-assembly and I know the MG reads look sound. Unlike the HMP which only had like 10 species per sample for some reason


# Create table of gene callers + assigned taxonomy
cut -f 1,3 /scratch/midway3/sjarif/07_MT_exports_HMP/MSM79H54_kaiju.out > /scratch/midway3/sjarif/07_MT_exports_HMP/MSM79H54_gene_taxonomy.tsv
# column 1 = gene ID, column 2 = species name or NCBI taxon ID

# Import taxonomy
anvi-import-taxonomy-for-genes \
    -i /scratch/midway3/sjarif/07_MT_exports_HMP/MSM79H54_gene_taxonomy.tsv \
    -c /scratch/midway3/sjarif/07_MT_exports_HMP/MSM79H54_contigs.db \
    --namespace Kaiju











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

# Total genes (Prodigal): 4,479,066
# Genes kept (≥50 aa AND non-partial):~1,456,735
# Genes removed: 3,022,331

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


# Run kaiju for MM
#!/bin/bash -e
#SBATCH -t 36:00:00 -N 1
#SBATCH --mem=150G 
#SBATCH --ntasks=1
#SBATCH -p blekhman
#SBATCH --account pi-blekhman
#SBATCH --mail-type=END
#SBATCH --mail-user=sjarif@uchicago.edu

kaiju -t /scratch/midway3/sjarif/db/kaijudb/nodes.dmp \
      -f /scratch/midway3/sjarif/db/kaijudb/progenomes/kaiju_db_progenomes.fmi \
      -i /scratch/midway3/sjarif/07_MT_exports_HMP/masiManimba_genes.faa \
      -o /scratch/midway3/sjarif/07_MT_exports_HMP/masiManimba_kaiju.out \
      -z 16 \
      -a greedy \
      -e 5 \
      -p \
      -v


# still zero hits
#!/bin/bash -e
#SBATCH -t 36:00:00 -N 1
#SBATCH --mem=300G 
#SBATCH --ntasks=1
#SBATCH -p bigmem
#SBATCH --account pi-blekhman
#SBATCH --mail-type=END
#SBATCH --mail-user=sjarif@uchicago.edu

cd /scratch/midway3/sjarif/db/kaiju_refseqnr
kaiju-makedb -s refseq_nr


# trying with diamond because kaiju is giving me issues
# (works!)
#!/bin/bash -e
#SBATCH -t 36:00:00 -N 1
#SBATCH --mem=50G 
#SBATCH --ntasks=1
#SBATCH -p blekhman
#SBATCH --account pi-blekhman
#SBATCH --mail-type=END
#SBATCH --mail-user=sjarif@uchicago.edu

diamond makedb \
  --in refseq_ref.faa \
  --db refseq_ref


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


# here is the output

head masiManimba_diamond_matches.tsv

4194306 WP_153206372.1  55.2    315     3.46e-102       313
4194306 WP_089290180.1  54.9    315     1.98e-101       311
4194306 WP_354080009.1  55.6    313     2.71e-101       310
4194306 WP_096702232.1  52.3    321     1.22e-98        303
4194306 WP_075771251.1  53.7    315     2.85e-98        302
4194306 WP_011384511.1  52.2    320     9.24e-98        301
4194306 WP_343865446.1  54.1    316     2.18e-97        300
4194306 WP_309420228.1  54.1    314     2.46e-97        300
4194306 WP_101248786.1  54.8    314     3.61e-97        300
4194306 WP_114394883.1  53.2    316     3.73e-97        300


# now filter for the best hit (lowest e-value first, highest bitscore first)
# this does not filter by 90% identity for example, but can do that later in R
# if we get some interesting gene expression patterns we can return to the large tsv
sort -k1,1 -k5,5g -k6,6nr masiManimba_diamond_matches.tsv \
  | awk '!seen[$1]++' \
  > masiManimba_diamond_matches_besthit.tsv

# finally, add taxonomy
# get list of accessions
cut -f2 masiManimba_diamond_matches_besthit.tsv | sort -u > masiManimba_besthit_accessions.txt

# map to accession IDs to taxonomy IDs
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz

zgrep -F -f masiManimba_besthit_accessions.txt ../db/diamond/prot.accession2taxid.gz \
  | awk '{print $1 "\t" $3}' \
  > masiManimba_accession_to_taxid.tsv

# download taxonomy
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
mkdir ncbi_taxdump
tar -xzf taxdump.tar.gz -C ncbi_taxdump

# extract parent + rank
awk -F '\t\\|\t' '{print $1"\t"$2"\t"$3}' \
  ../db/diamond/ncbi_taxdump/nodes.dmp \
  > ../db/diamond/ncbi_taxdump/taxid_nodes.tsv

# extract scientific names
awk -F '|' '$4 ~ /scientific name/ {print $1"\t"$2}' \
  ../db/diamond/ncbi_taxdump/names.dmp \
  > ../db/diamond/ncbi_taxdump/taxid_names.tsv


# check structure
head /scratch/midway3/sjarif/db/diamond/ncbi_taxdump/taxid_nodes.tsv
1       1       no rank
2       131567  domain
6       335928  genus
7       6       species
9       32199   species
10      1706371 genus
11      1707    species
13      203488  genus
14      13      species
16      32011   genus

head /scratch/midway3/sjarif/db/diamond/ncbi_taxdump/taxid_names.tsv
1                       root
2                       Bacteria
6                       Azorhizobium
7                       Azorhizobium caulinodans
9                       Buchnera aphidicola
10                      Cellvibrio
11                      Cellulomonas gilvus
13                      Dictyoglomus
14                      Dictyoglomus thermophilum
16                      Methylophilus



# python script to extract LONG lineage
# (works!!)
#!/usr/bin/env python3

nodes_file = "taxid_nodes.tsv"
names_file = "taxid_names.tsv"
output_file = "taxid_lineage_long.tsv"


names = {}
with open(names_file) as f:
    for line in f:
        parts = line.rstrip().split("\t", 1)  # split only on first tab
        if len(parts) == 2:
            taxid, name = parts
            names[taxid.strip()] = name.strip()
        else:
            continue  # skip malformed lines


with open(output_file, "w") as out:
    out.write("taxid\trank\tname\n")
    with open(nodes_file) as f:
        for line in f:
            parts = line.rstrip().split("\t")
            if len(parts) == 3:
                taxid, parent, rank = parts
                name = names.get(taxid.strip(), "NA")
                out.write(f"{taxid.strip()}\t{rank.strip()}\t{name}\n")

# check format

# TAXID + RANK INFORMATION AND END-OF-NODE NAME
head /scratch/midway3/sjarif/db/diamond/ncbi_taxdump/taxid_lineage_long.tsv
taxid   rank    name
1       no rank root
2       domain  Bacteria
6       genus   Azorhizobium
7       species Azorhizobium caulinodans
9       species Buchnera aphidicola
10      genus   Cellvibrio
11      species Cellulomonas gilvus
13      genus   Dictyoglomus
14      species Dictyoglomus thermophilum

# TAXID + NCBI ACCESSIONS
head /scratch/midway3/sjarif/07_MT_exports_HMP/masiManimba_accession_to_taxid.tsv
WP_000002019    1301
WP_000002117    543
WP_000002283    543
WP_000002446    543
WP_000002542    543
WP_000002800    2
WP_000002830    562
WP_000002907    543
WP_000002953    1236
WP_000003060    543

# MASI MANIMBA GENE CALLS, CONTIG IDS AND AA SEQUENCES
head /scratch/midway3/sjarif/07_MT_exports_HMP/masiManimba_gene_calls.tsv
gene_callers_id contig  start   stop    direction       partial call_type       source  version aa_sequence
0       c_000000000001  2       335     f       1       1       prodigal        v2.6.3  VMLMTFFVEHYGYYFDQELLAQSAGYVRDALVMACLDNYSEYEYLERILQDAICTEPIAETFAEERPASISEKYQKYQSKHYEPAPHEYVEYKTKNTYSEDPIAGKASKK
1       c_000000000009  0       612     r       1       1       prodigal        v2.6.3  TEREFIPVLLGGDINAYSVARAFYEEYQVKSLVFGKYQTGPAYRSQIIDYTPNVDIDTMPVMLKTVNDIARAHADKTIVLMGCGDNYVALVAQAKDANELADNIVAPYAPYSMLEQCQKKEIFYELCEKHGVPYPHTFTFTAQMLNAQGEAPAEVLDQIDFPFPMILKPSDGIMWWQHEFEGQKKAYEIADRAELEQVIHDSYA
2       c_000000000017  26      551     f       0       1       prodigal        v2.6.3  MAKYRKLGRTSSQRKALIRSQVTALLHNGRIITTEARAQEVRKVAEGLIASAVKECDNFEEVTVKAKVARKDSEGKRVKEVVDGKKVTVYDEVEKTIKKDMPSRLHARREMLKVLYPVTETGAKKKDTKEVDLVDKLFTEIAPKYKDRNGGYTRIVKIGLRKGDGAMEVVLELV
3       c_000000000025  46      664     f       1       1       prodigal        v2.6.3  MEKKLQTLKNILDSSRYTVALCGSGILKECGYPGILSLDVAYDIENRYGDSPEYIYSSAYFSTRPEKFFKFYKSEILDKDLEPSETFYALAELEKQEKLQTIISKNSFSLSERAGCKHVYNIYGTIHKNICTHCGKEYPVEFVKTSPSVPLCEECGHIIRPNVKLFGEMLDHEIIAKIADEIAKADVLLLLGTSLESYTYKNYIKY

# MASI MANIMBA GENE CALLS AND AA SEQUENCES
head /scratch/midway3/sjarif/07_MT_exports_HMP/masiManimba_genes.faa
>4194306
MQNLGLIADIGGTNARFALCAADGTVERPAVLACRDYSGLADAVSAYCKQAGVAPETVKTAAFAVACPVLGDKVSLTNCAWAFSIEELKTKLGLDNLTVVNDFVAQALAVPELKESEKVK
IGRGEPQAGHPIAVVGPGTGLGVSILVPQNDGAWTALPSEGGHATMPAPYPEEERVISALRGEYGHVSAERVVSGMGLTNLYKTLKVLRNEPADVLPPEEITNRALAGDALCREAVQMMF
GLLGTLAGNLALTVGALGGVYIAGGIVPREGLLEMFKESAFRVRFEAKGRFTDYLSRIPTYVMTAEYPAFLGLSALIKKNFKK
>4194307
MMWDDFLPLGQALEARFPETDVLSVSDAELKRMLAALPMADGAEAMPESADFFYKIKVAWIQARRTEPADTSMEADL
>2
MAKYRKLGRTSSQRKALIRSQVTALLHNGRIITTEARAQEVRKVAEGLIASAVKECDNFEEVTVKAKVARKDSEGKRVKEVVDGKKVTVYDEVEKTIKKDMPSRLHARREMLKVLYPVTE
TGAKKKDTKEVDLVDKLFTEIAPKYKDRNGGYTRIVKIGLRKGDGAMEVVLELV
>4194312







MasimManimba contigs (scg.db)
        │
        ▼
Gene calls (masiManimba_gene_calls.tsv)
        │
        ▼
Filtered gene IDs (masiManimba_gene_ids_min50aa_complete.txt)
        │
        ▼
AA sequences (masiManimba_genes.faa)
        │
        ▼
DIAMOND search → Best hits (masiManimba_diamond_matches_besthit.tsv)
        │
        ▼
Map DIAMOND protein ID → NCBI taxid (masiManimba_accession_to_taxid.tsv)
        │
        ▼
Attach taxonomy lineage (taxid_lineage_long.tsv)



# ANVIO GENE CALLERS AS ROWS, MTX SAMPLES AS COLUMNS: gene coverage
/project/blekhman/sjarif/konzo/07_MT_exports/gene_coverages/merged_gene_coverages.tsv

# ANVIO GENE CALLERS AS ROWS, MTX SAMPLES AS COLUMNS: gene detection
/project/blekhman/sjarif/konzo/07_MT_exports/gene_detection/merged_gene_detection.tsv

# ANVIO GENE CALLERS AS ROWS X NCBI PROTEIN ACCESSION (one per gene; best hit)
/scratch/midway3/sjarif/07_MT_exports_HMP/masiManimba_diamond_matches_besthit.tsv

# NCBI PROTEIN ACCESSION X NCBI TAXID
/scratch/midway3/sjarif/07_MT_exports_HMP/masiManimba_accession_to_taxid.tsv

# NCBI TAXID + RANK INFORMATION AND END-OF-NODE NAME
/scratch/midway3/sjarif/db/diamond/ncbi_taxdump/taxid_lineage_long.tsv

