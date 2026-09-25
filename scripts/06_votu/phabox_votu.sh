#!/bin/bash
#SBATCH --job-name=phabox_votu
#SBATCH --output=logs/phabox_votu_%j.out
#SBATCH --error=logs/phabox_votu_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=batch
#SBATCH --account=CHANGE_ME
#
# Species-level vOTU clustering (95% ANI / 85% coverage — the community-
# standard threshold used by most published viral catalogs) on the
# BLAST-safe main catalog fasta (05_catalog_build/fix_headers.py output).
# This is a single all-vs-all run, not a per-sample array job.

source activate phabox_env
source "$(dirname "$0")/../../config.sh"

phabox2 --task votu \
    --dbdir "${PHABOX_DB}" \
    --outpth "${VOTU_OUT}" \
    --contigs "${CATALOG_MAIN}_safe.fasta" \
    --threads 16 \
    --mode ANI \
    --ani 95 \
    --tcov 85
