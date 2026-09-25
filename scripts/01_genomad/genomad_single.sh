#!/bin/bash
#SBATCH --job-name=genomad_single
#SBATCH --output=logs/genomad_single_%j.out
#SBATCH --error=logs/genomad_single_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=batch
#SBATCH --account=CHANGE_ME
#
# Reruns geNomad on ONE named contig file. Use this for reprocessing
# specific failed/missing samples instead of resubmitting an array job
# with guessed index numbers — indices can silently point at the wrong
# file if the contigs directory changed since the array was built.
#
# Usage: sbatch genomad_single.sh /path/to/contigs/sample.fa

module load genomad/1.7.0
source "$(dirname "$0")/../../config.sh"

FASTA_FILE="$1"
if [ -z "${FASTA_FILE}" ]; then
    echo "Usage: sbatch genomad_single.sh <path-to-fasta>" >&2
    exit 1
fi

FILE_NAME=$(basename "${FASTA_FILE}")
SAMPLE_NAME="${FILE_NAME%.*}"
SAMPLE_OUT_DIR="${GENOMAD_OUT}/${SAMPLE_NAME}"
mkdir -p "${SAMPLE_OUT_DIR}"

genomad end-to-end --threads 16 "${FASTA_FILE}" "${SAMPLE_OUT_DIR}" "${GENOMAD_DB}"
