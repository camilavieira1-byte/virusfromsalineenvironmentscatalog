#!/bin/bash
#SBATCH --job-name=phabox_contam_genomad
#SBATCH --output=logs/phabox_contam_genomad_%A_%a.out
#SBATCH --error=logs/phabox_contam_genomad_%A_%a.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=batch
#SBATCH --account=CHANGE_ME
# sbatch --array=1-$(wc -l < virus_fasta_list.txt)%10 phabox_contamination_from_genomad.sh
#
# Runs phabox2's contamination/provirus task on geNomad's identified
# viral sequences. Writes to its OWN output root (separate from
# phabox_contamination_from_phabox.sh) so the two checks never overwrite
# each other's final_prediction/contamination_prediction.tsv.

source activate phabox_env
source "$(dirname "$0")/../../config.sh"

FASTA_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${GENOMAD_VIRUS_LIST}")
SAMPLE_NAME=$(basename "${FASTA_FILE}" _virus.fna)
SAMPLE_OUT_DIR="${PHABOX_CONTAM_GENOMAD_OUT}/${SAMPLE_NAME}"
mkdir -p "${SAMPLE_OUT_DIR}"

phabox2 --task contamination \
    --dbdir "${PHABOX_DB}" \
    --outpth "${SAMPLE_OUT_DIR}" \
    --contigs "${FASTA_FILE}" \
    --threads 16
