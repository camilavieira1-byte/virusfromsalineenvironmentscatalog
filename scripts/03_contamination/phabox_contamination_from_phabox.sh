#!/bin/bash
#SBATCH --job-name=phabox_contam_phabox
#SBATCH --output=logs/phabox_contam_phabox_%A_%a.out
#SBATCH --error=logs/phabox_contam_phabox_%A_%a.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=batch
#SBATCH --account=CHANGE_ME
# sbatch --array=1-$(wc -l < phabox_uncertain_list.txt)%10 phabox_contamination_from_phabox.sh
#
# Runs the contamination task on phabox2's OWN uncertain-sequence bucket
# (its documented, intended input for this task — see phabox2's wiki).
# This extends the SAME per-sample output folder phabox2's end_to_end run
# already created (adds contamination_prediction.tsv alongside the
# existing phamer/phagcn/phatyp/cherry results), rather than making a
# new output tree.

source activate phabox_env
source "$(dirname "$0")/../../config.sh"

FASTA_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${PHABOX_UNCERTAIN_LIST}")
# FASTA_FILE looks like: .../phabox_outputs/<sample>/final_prediction/phamer_supplementary/uncertain_sequences_for_contamination_task.fa
SAMPLE_OUT_DIR=$(dirname "$(dirname "$(dirname "${FASTA_FILE}")")")

phabox2 --task contamination \
    --dbdir "${PHABOX_DB}" \
    --outpth "${SAMPLE_OUT_DIR}" \
    --contigs "${FASTA_FILE}" \
    --threads 16
