#!/bin/bash
#SBATCH --job-name=mergetrim_main
#SBATCH --output=logs/mergetrim_main_%j.out
#SBATCH --error=logs/mergetrim_main_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --partition=batch
#SBATCH --account=CHANGE_ME

source activate phabox_env   # or any env with pandas installed
source "$(dirname "$0")/../../config.sh"

python "$(dirname "$0")/build_comparison_table.py"
python "$(dirname "$0")/mergetrim_main.py"
python "$(dirname "$0")/fix_headers.py" "${CATALOG_MAIN}" MAIN
