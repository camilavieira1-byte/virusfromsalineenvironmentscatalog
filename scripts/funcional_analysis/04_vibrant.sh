#!/bin/bash
#SBATCH --job-name=vibrant_amg_complete
#SBATCH --output=logs/vibrant_amg_complete_%j.out
#SBATCH --error=logs/vibrant_amg_complete_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=batch
##SBATCH --account=<your-account>   # uncomment and edit if your cluster requires it

# Run VIBRANT on the vOTU set. Run from the project directory, after `mkdir -p logs`.
#
# Configurable via environment variables:
#   PROJECT_DIR  working directory (default: directory where sbatch was called)
#   INPUT_FASTA  vOTU fasta        (default: $PROJECT_DIR/complete_vOTUs.fasta)
#   OUT_DIR      output directory  (default: $PROJECT_DIR/vibrant_amg_complete_output)

set -eo pipefail

PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
INPUT_FASTA="${INPUT_FASTA:-$PROJECT_DIR/complete_vOTUs.fasta}"
OUT_DIR="${OUT_DIR:-$PROJECT_DIR/vibrant_amg_complete_output}"

module load vibrant/1.2.1   # adjust for your cluster
cd "$PROJECT_DIR"

VIBRANT_run.py -i "$INPUT_FASTA" -folder "$OUT_DIR" -t "${SLURM_CPUS_PER_TASK:-16}"
