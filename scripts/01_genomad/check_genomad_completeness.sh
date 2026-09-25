#!/bin/bash
# Prints which contig samples do NOT yet have a completed geNomad run
# (checks for the actual summary output, not just an output folder —
# a folder can exist from a crashed/incomplete run).
set -euo pipefail
source "$(dirname "$0")/../../config.sh"

comm -23 \
  <(ls "${CONTIGS_DIR}"/*.fa "${CONTIGS_DIR}"/*.fasta 2>/dev/null | xargs -n1 basename | sed -E 's/\.(fa|fasta)$//' | sort -u) \
  <(find "${GENOMAD_OUT}" -name "*_virus_summary.tsv" | xargs -n1 dirname | xargs -n1 dirname | xargs -n1 basename | sort -u)
