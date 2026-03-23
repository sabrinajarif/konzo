**Recommended steps for next user as of March 23, 2026**
- Here, I used only a small subset of the Democratic Republic of Congo samples: the Masi Manimba samples. For more power, expand to other regions as well.
- I used Human Microbiome Project samples as an "industrialized" comparison. There were several issues with this dataset (i. participants were much older, ii. assembly was performed at the single-sample level which meant iii. only a few high abundant species were detected per sample), and the next user should either pick a better contrasting sample set or just not use HMP.
- Metagenomic assembly was performed in 2023/2024, and steps should probably be reviewed and updated. Co-assembly is recommended, but if too computationaly expensive, limit co-assembly to each geographic region. 
- Every step in the preprocessing (01_, 02_, 03_ etc) contains its own Snakemake file. This was to help me troubleshoot- ideally this should be streamlined into one or two Snaekmake workflows. 
- I used QC'd metatranscriptomic and metagenomic samples to begin this workflow. I noticed that the QC'd MT reads that were available were still heavily contaminated with rRNA after a BLAST search. The next user should update and optimize the QC starting with the raw reads; the only additional QC I performed was utilizing Ribodetector for the MTs. 
- Ideally, the metagemonic contig library should be annotated with taxonomic information at the genome level as well as at the gene level. There should also be functional annotation (KEGG, etc) at the gene level. Because I needed this data ASAP, I just used DIAMOND to annotate genes at the taxonomic level in the MT data after aligning to contigs, choosing the best hit for each. The next user should try different approaches for annotation, and ensure the MG and MT data are directly comparable.
- Since I just took ratios of RNA:DNA for my analysis, I did not transform the counts in any way. The next user should optimize this/ check what is current in MT literature.

**Directory structure**

*metadata*\
-> 

*preprocessing*\
->

*analysis*\
->
