#!/bin/bash
#SBATCH --job-name=virsorter2_dramv_complete
#SBATCH --output=logs/virsorter2_dramv_complete_%j.out
#SBATCH --error=logs/virsorter2_dramv_complete_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=batch
##SBATCH --account=<your-account>   # uncomment and edit if your cluster requires it

# Run VirSorter2 on the vOTU set and prepare the files DRAM-v needs.
# Run from the project directory, after `mkdir -p logs`.
#
# Configurable via environment variables:
#   PROJECT_DIR  working directory            (default: directory where sbatch was called)
#   INPUT_FASTA  vOTU fasta                   (default: $PROJECT_DIR/complete_vOTUs.fasta)
#   VS2_SIF      VirSorter2 Singularity image (default: $PROJECT_DIR/virsorter2.sif)
#   VS2_DB       VirSorter2 database          (default: $PROJECT_DIR/virsorter2_db)
#   OUT_DIR      output directory             (default: $PROJECT_DIR/virsorter2_for_dramv_complete)
#   BIND_PATHS   extra Singularity bind paths, comma-separated (e.g. a database on another filesystem)

set -eo pipefail

PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
INPUT_FASTA="${INPUT_FASTA:-$PROJECT_DIR/complete_vOTUs.fasta}"
VS2_SIF="${VS2_SIF:-$PROJECT_DIR/virsorter2.sif}"
VS2_DB="${VS2_DB:-$PROJECT_DIR/virsorter2_db}"
OUT_DIR="${OUT_DIR:-$PROJECT_DIR/virsorter2_for_dramv_complete}"
BIND="$PROJECT_DIR:$PROJECT_DIR:rw${BIND_PATHS:+,$BIND_PATHS}"

module load singularity   # adjust for your cluster

singularity exec \
    --bind "$BIND" \
    "$VS2_SIF" \
    virsorter run \
    -i "$INPUT_FASTA" \
    -w "$OUT_DIR" \
    -d "$VS2_DB" \
    --prep-for-dramv \
    -j "${SLURM_CPUS_PER_TASK:-16}" \
    all
