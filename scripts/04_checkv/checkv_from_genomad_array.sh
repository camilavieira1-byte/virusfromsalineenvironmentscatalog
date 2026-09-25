#!/bin/bash
#SBATCH --job-name=checkv_genomad
#SBATCH --output=logs/checkv_genomad_%A_%a.out
#SBATCH --error=logs/checkv_genomad_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=batch
#SBATCH --account=CHANGE_ME
# sbatch --array=1-$(wc -l < virus_fasta_list.txt)%10 checkv_from_genomad_array.sh

module load checkv/1.0.1
source "$(dirname "$0")/../../config.sh"

# IMPORTANT: some HPC checkv modules export the database path under a
# DIFFERENT env var name than checkv itself expects (e.g. $CHECKV_DB
# instead of $CHECKVDB). Check with `echo $CHECKVDB` right after your
# `module load` line; if it's empty, uncomment and set this explicitly:
# export CHECKVDB="/path/to/checkv-db-vX.X"

FASTA_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${GENOMAD_VIRUS_LIST}")
SAMPLE_NAME=$(basename "${FASTA_FILE}" _virus.fna)
SAMPLE_OUT_DIR="${CHECKV_GENOMAD_OUT}/${SAMPLE_NAME}"
mkdir -p "${SAMPLE_OUT_DIR}"

checkv end_to_end "${FASTA_FILE}" "${SAMPLE_OUT_DIR}" -t 16
