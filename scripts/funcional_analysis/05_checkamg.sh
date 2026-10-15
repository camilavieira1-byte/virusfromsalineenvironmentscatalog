#!/bin/bash
#SBATCH --job-name=checkamg_complete
#SBATCH --output=logs/checkamg_complete_%j.out
#SBATCH --error=logs/checkamg_complete_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=30
#SBATCH --mem=200G
#SBATCH --partition=batch
##SBATCH --account=<your-account>   # uncomment and edit if your cluster requires it

# Run CheckAMG end-to-end on the vOTU set. Run from the project directory, after `mkdir -p logs`.
#
# Configurable via environment variables:
#   CHECKAMG_DB  directory holding the CheckAMG annotate/de-novo databases (REQUIRED)
#   PROJECT_DIR  working directory (default: directory where sbatch was called)
#   INPUT_FASTA  vOTU fasta        (default: $PROJECT_DIR/complete_vOTUs.fasta)
#   OUT_DIR      output directory  (default: $PROJECT_DIR/checkamg_complete_out)
#   CONDA_SH     path to conda.sh  (default: $HOME/miniconda3/etc/profile.d/conda.sh)
#   CONDA_ENV    conda env with CheckAMG installed (default: checkamg_env)
#   MEM_GB       memory passed to CheckAMG with -m (default: 180, keep below --mem)

set -eo pipefail

: "${CHECKAMG_DB:?Set CHECKAMG_DB to the CheckAMG database directory}"
PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
INPUT_FASTA="${INPUT_FASTA:-$PROJECT_DIR/complete_vOTUs.fasta}"
OUT_DIR="${OUT_DIR:-$PROJECT_DIR/checkamg_complete_out}"
CONDA_SH="${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-checkamg_env}"
MEM_GB="${MEM_GB:-180}"

# shellcheck disable=SC1090
source "$CONDA_SH"
conda activate "$CONDA_ENV"
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

cd "$PROJECT_DIR"
mkdir -p "$OUT_DIR"

checkamg end-to-end \
  -d "$CHECKAMG_DB" \
  -i "$INPUT_FASTA" \
  -o "$OUT_DIR" \
  --input-type nucl \
  --denovo-accelerator cpu \
  -t "$SLURM_CPUS_PER_TASK" \
  -m "$MEM_GB"
