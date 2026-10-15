#!/bin/bash
#SBATCH --job-name=dramv_distill
#SBATCH --output=dramv_distill_%j.out
#SBATCH --error=dramv_distill_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=batch
##SBATCH --account=<your-account>   # uncomment and edit if your cluster requires it

# Distill DRAM-v annotations (run 02_dramv_annotate.sh first).
#
# Configurable via environment variables:
#   PROJECT_DIR   working directory (default: directory where sbatch was called)
#   ANNOTATIONS   annotations table to distill
#                 (default: $PROJECT_DIR/dramv_annotate_out/annotations_with_amg_flags.tsv)
#   OUT_DIR       distill output (default: $PROJECT_DIR/dramv_distill_out)

set -eo pipefail

PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
ANNOTATIONS="${ANNOTATIONS:-$PROJECT_DIR/dramv_annotate_out/annotations_with_amg_flags.tsv}"
OUT_DIR="${OUT_DIR:-$PROJECT_DIR/dramv_distill_out}"

module load dram/1.5.0   # adjust for your cluster

DRAM-v.py distill \
  -i "$ANNOTATIONS" \
  -o "$OUT_DIR"
