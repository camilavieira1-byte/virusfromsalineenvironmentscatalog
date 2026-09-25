---

### `config.sh`
```bash
#!/bin/bash
# ============================================================
# Central configuration for the viral catalog pipeline.
# Edit the paths below for your own environment, then source
# this file near the top of every script, e.g.:
#   source "$(dirname "$0")/../../config.sh"
# (adjust the relative path depth to match where the script lives)
# ============================================================

# Root directory for the whole project
export BASE_DIR="/path/to/your/project"

# Input contigs: one .fa/.fasta file per sample or coassembly
export CONTIGS_DIR="${BASE_DIR}/contigs"

# Reference databases — download separately per each tool's own instructions
export GENOMAD_DB="${BASE_DIR}/genomad_db"
export PHABOX_DB="${BASE_DIR}/phabox_db_v2_2"
# checkV's db path is often exported by its module under a DIFFERENT name
# than checkv itself expects (e.g. $CHECKV_DB instead of$CHECKVDB) — the
# 04_checkv scripts have a commented-out line to set this explicitly if needed.

# Identification step outputs (Step 01)
export GENOMAD_OUT="${BASE_DIR}/genomad_outputs"
export PHABOX_OUT="${BASE_DIR}/phabox_outputs"                     # end_to_end output (phamer/phagcn/phatyp/cherry/phavip)

# --- Static sample manifests ---
# SLURM array scripts read a fixed line from these files via
# `sed -n "${SLURM_ARRAY_TASK_ID}p"` rather than indexing into a live
# `contigs/*.fa` glob. This matters: a live glob's order can shift if
# files are added/removed between when you compute an index and when
# the array job actually runs, silently processing the wrong sample.
# Regenerate manifests (scripts below) and resubmit by array size, not
# by reusing old index numbers, whenever contigs/ or an output dir changes.
export CONTIG_MANIFEST="${BASE_DIR}/contig_manifest.txt"
export GENOMAD_VIRUS_LIST="${BASE_DIR}/virus_fasta_list.txt"
export PHABOX_PREDICTED_VIRUS_LIST="${BASE_DIR}/phabox_predicted_virus_list.txt"
export PHABOX_UNCERTAIN_LIST="${BASE_DIR}/phabox_uncertain_list.txt"

# Contamination check outputs (Step 03)
export PHABOX_CONTAM_GENOMAD_OUT="${BASE_DIR}/phabox_contamination_genomad"

# checkV outputs (Step 04)
export CHECKV_GENOMAD_OUT="${BASE_DIR}/checkv_genomad_outputs"
export CHECKV_PHABOX_OUT="${BASE_DIR}/checkv_phabox_outputs"

# Non-redundant catalogs (basenames — .fasta/.tsv appended by the scripts)
export CATALOG_MAIN="${BASE_DIR}/catalog_main_nonredundant"
export CATALOG_RELAXED="${BASE_DIR}/catalog_relaxed_nonredundant"

# vOTU clustering + annotation (Step 06)
export VOTU_OUT="${BASE_DIR}/votu_main_catalog"
export VOTU_ANNOTATED_OUT="${BASE_DIR}/votu_catalog_annotated"

mkdir -p "${BASE_DIR}/logs"
