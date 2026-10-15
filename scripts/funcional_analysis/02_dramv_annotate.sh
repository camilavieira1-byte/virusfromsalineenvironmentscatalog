#!/bin/bash
#SBATCH --job-name=dramv_complete
#SBATCH --output=dramv_complete_%j.out
#SBATCH --error=dramv_complete_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --partition=batch
##SBATCH --account=<your-account>   # uncomment and edit if your cluster requires it

# Annotate the VirSorter2 output with DRAM-v (run 01_virsorter2_for_dramv.sh first).
#
# Configurable via environment variables:
#   PROJECT_DIR  working directory (default: directory where sbatch was called)
#   VS2_DIR      VirSorter2 output (default: $PROJECT_DIR/virsorter2_for_dramv_complete)
#   OUT_DIR      annotation output (default: $PROJECT_DIR/dramv_annotate_complete_out)
#   MIN_CONTIG   minimum contig size in bp (default: 2500)

set -eo pipefail

PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
VS2_DIR="${VS2_DIR:-$PROJECT_DIR/virsorter2_for_dramv_complete}"
OUT_DIR="${OUT_DIR:-$PROJECT_DIR/dramv_annotate_complete_out}"
MIN_CONTIG="${MIN_CONTIG:-2500}"

module load dram/1.5.0   # adjust for your cluster
cd "$PROJECT_DIR"

DRAM-v.py annotate \
  -i "$VS2_DIR/for-dramv/final-viral-combined-for-dramv.fa" \
  -v "$VS2_DIR/for-dramv/viral-affi-contigs-for-dramv.tab" \
  -o "$OUT_DIR" \
  --min_contig_size "$MIN_CONTIG" \
  --threads "${SLURM_CPUS_PER_TASK:-32}"
