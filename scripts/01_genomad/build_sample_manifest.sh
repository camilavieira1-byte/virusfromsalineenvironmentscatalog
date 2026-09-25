#!/bin/bash
# Generates a stable, one-per-line manifest of contig fasta files.
# Run this once before submitting the geNomad array job. If contigs/
# gains or loses files later, regenerate this file and resubmit any
# remaining/failed samples BY NAME (genomad_single.sh), not by reusing
# old array index numbers — the whole point of the manifest is that
# indices only stay meaningful for the manifest that produced them.
set -euo pipefail
source "$(dirname "$0")/../../config.sh"

find "${CONTIGS_DIR}" -maxdepth 1 \( -name "*.fa" -o -name "*.fasta" \) | sort > "${CONTIG_MANIFEST}"
echo "Wrote $(wc -l < "${CONTIG_MANIFEST}") entries to ${CONTIG_MANIFEST}"
