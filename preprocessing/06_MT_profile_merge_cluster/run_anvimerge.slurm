#!/bin/bash
#SBATCH --job-name=anvi_merge
#SBATCH --output=anvi_merge.out
#SBATCH --error=anvi_merge.err
#SBATCH --time=12:00:00
#SBATCH --partition=blekhman
#SBATCH --account=pi-blekhman
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

cd /project/blekhman/sjarif/konzo
conda activate anvio-dev

BASE="/project/blekhman/sjarif/konzo/06_MT_profile_merge_cluster"
MT_PROFILES="$BASE/profiled/mtx"
MG_PROFILES="$BASE/profiled/mgx"
CONTIGS="/project/blekhman/shared/konzo/contigs/masiManimba.db"

# ----- MT PROFILE.DBS -----
MT_FILES=($MT_PROFILES/*/PROFILE.db)

# ----- MG PROFILE.DBS -----
MG_FILES=($MG_PROFILES/*/PROFILE.db)

### ----- RUN ANVI-MERGE ----- ###
anvi-merge \
    $MT_FILES \
    $MG_FILES \
    -c "$CONTIGS" \
    -o "$BASE/merged2" \
    --skip-hierarchical-clustering

echo "Done."
