#!/bin/bash
#SBATCH --job-name=phabox_annotate
#SBATCH --output=logs/phabox_annotate_%j.out
#SBATCH --error=logs/phabox_annotate_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=batch
#SBATCH --account=CHANGE_ME
#
# Taxonomy (phagcn), lifestyle (phatyp), host prediction (cherry), and
# protein annotation (phavip) on the FINAL vOTU representative sequences
# (one per species-level cluster), not the raw pre-clustering catalog.
#
# --skip Y is important: without it, phabox2 re-runs virus identification
# (phamer) on sequences already confirmed viral by this whole pipeline,
# which is redundant and can even mislabel a legitimate but divergent
# virus as low-confidence.

source activate phabox_env
source "$(dirname "$0")/../../config.sh"

phabox2 --task end_to_end \
    --dbdir "${PHABOX_DB}" \
    --outpth "${VOTU_ANNOTATED_OUT}" \
    --contigs "${VOTU_OUT}/final_prediction/votu_supplementary/ANI_based_vOTU.fasta" \
    --threads 16 \
    --skip Y
