#!/bin/bash
# Builds the static per-sample fasta lists consumed by every downstream
# array job (contamination task, checkV). Regenerate this any time new
# geNomad or phabox2 runs complete, BEFORE resubmitting downstream steps,
# so array indices stay aligned with what's actually on disk.
set -euo pipefail
source "$(dirname "$0")/../../config.sh"

# geNomad-identified viral sequences (skip samples where geNomad found none)
find "${GENOMAD_OUT}" -name "*_virus.fna" -size +0c | sort > "${GENOMAD_VIRUS_LIST}"
echo "geNomad virus fastas: $(wc -l < "${GENOMAD_VIRUS_LIST}")"

# phabox2's own confidently-identified viral sequences (from its phamer step)
find "${PHABOX_OUT}" -name "predicted_virus.fa" -size +0c | sort > "${PHABOX_PREDICTED_VIRUS_LIST}"
echo "phabox predicted_virus fastas: $(wc -l < "${PHABOX_PREDICTED_VIRUS_LIST}")"

# phabox2's own "uncertain" bucket — the file phabox2 itself designates as
# the intended input for its contamination task (NOT predicted_virus.fa)
find "${PHABOX_OUT}" -name "uncertain_sequences_for_contamination_task.fa" -size +0c | sort > "${PHABOX_UNCERTAIN_LIST}"
echo "phabox uncertain-sequence fastas: $(wc -l < "${PHABOX_UNCERTAIN_LIST}")"
