#!/bin/bash
#SBATCH --job-name=checkv_phabox
#SBATCH --output=logs/checkv_phabox_%A_%a.out
#SBATCH --error=logs/checkv_phabox_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=batch
#SBATCH --account=CHANGE_ME
# sbatch --array=1-$(wc -l < phabox_predicted_virus_list.txt)%10 checkv_from_phabox_array.sh

module load checkv/1.0.1
source "$(dirname "$0")/../../config.sh"
# export CHECKVDB="/path/to/checkv-db-vX.X"   # see note in checkv_from_genomad_array.sh

FASTA_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${PHABOX_PREDICTED_VIRUS_LIST}")
# FASTA_FILE looks like: .../phabox_outputs/<sample>/final_prediction/phamer_supplementary/predicted_virus.fa
SAMPLE_NAME=$(basename "$(dirname "$(dirname "$(dirname "${FASTA_FILE}")")")")
SAMPLE_OUT_DIR="${CHECKV_PHABOX_OUT}/${SAMPLE_NAME}"
mkdir -p "${SAMPLE_OUT_DIR}"

checkv end_to_end "${FASTA_FILE}" "${SAMPLE_OUT_DIR}" -t 16
