#!/bin/bash
#SBATCH --job-name=mergetrim_relaxed
#SBATCH --output=logs/mergetrim_relaxed_%j.out
#SBATCH --error=logs/mergetrim_relaxed_%j.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --partition=batch
#SBATCH --account=CHANGE_ME
# The relaxed catalog is ~100x larger than the main one (every quality
# tier, not just Complete/High/Medium) — give it real time and memory.

source activate phabox_env
source "$(dirname "$0")/../../config.sh"

# Requires comparison_genomad_pipeline_annotated.tsv — run
# mergetrim_main_slurm.sh (or at least build_comparison_table.py) first
# if that file doesn't exist yet.
python "$(dirname "$0")/mergetrim_relaxed.py"
python "$(dirname "$0")/fix_headers.py" "${CATALOG_RELAXED}" RELAXED
