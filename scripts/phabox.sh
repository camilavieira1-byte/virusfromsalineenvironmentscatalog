#!/bin/bash
#SBATCH --job-name=phabox2_array
#SBATCH --output=logs/phabox_%A_%a.out
#SBATCH --error=logs/phabox_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=batch
#SBATCH --account=CHANGE_ME
# Submit sized to the manifest, e.g.:
#   sbatch --array=1-$(wc -l < contig_manifest.txt)%10 phabox2_array.sh

# Load environment (uncomment your cluster's loading mechanism)
# module load phabox/2.0.0
# source activate phabox_env

source "$(dirname "$0")/../../config.sh"

FASTA_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${CONTIG_MANIFEST}")
if [ -z "${FASTA_FILE}" ] || [ ! -f "${FASTA_FILE}" ]; then
    echo "Error: no contig file for array index ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

FILE_NAME=$(basename "${FASTA_FILE}")
SAMPLE_NAME="${FILE_NAME%.*}"
SAMPLE_OUT_DIR="${PHABOX_OUT}/${SAMPLE_NAME}"
mkdir -p "${SAMPLE_OUT_DIR}"

echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Processing: ${FILE_NAME}"
echo "Output dir: ${SAMPLE_OUT_DIR}"

phabox2 \
    --task end_to_end \
    --dbdir "${PHABOX_DB}" \
    --contigs "${FASTA_FILE}" \
    --outpth "${SAMPLE_OUT_DIR}" \
    --threads "${SLURM_CPUS_PER_TASK:-16}" \
    --len 3000
