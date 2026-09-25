#!/bin/bash
#SBATCH --job-name=genomad_array
#SBATCH --output=logs/genomad_%A_%a.out
#SBATCH --error=logs/genomad_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=batch
#SBATCH --account=CHANGE_ME
# Submit sized to the manifest, e.g.:
#   sbatch --array=1-$(wc -l < contig_manifest.txt)%10 genomad_array.sh

module load genomad/1.7.0
source "$(dirname "$0")/../../config.sh"

FASTA_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${CONTIG_MANIFEST}")
if [ -z "${FASTA_FILE}" ] || [ ! -f "${FASTA_FILE}" ]; then
    echo "Error: no contig file for array index ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

FILE_NAME=$(basename "${FASTA_FILE}")
SAMPLE_NAME="${FILE_NAME%.*}"
SAMPLE_OUT_DIR="${GENOMAD_OUT}/${SAMPLE_NAME}"
mkdir -p "${SAMPLE_OUT_DIR}"

echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Processing: ${FILE_NAME}"
echo "Output dir: ${SAMPLE_OUT_DIR}"

genomad end-to-end \
    --threads 16 \
    "${FASTA_FILE}" \
    "${SAMPLE_OUT_DIR}" \
    "${GENOMAD_DB}"
