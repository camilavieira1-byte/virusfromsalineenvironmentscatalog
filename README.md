# Viral Catalog Pipeline

Builds a non-redundant, quality-controlled, species-level (vOTU) viral
genome catalog from metagenomic coassembly contigs, using two
independent identification tools (geNomad, phabox2) cross-validated
against each other and against checkV.

**Scope:** this covers identification → contamination/QC →
deduplication → vOTU clustering → taxonomy/host/lifestyle annotation.
It stops at a finished, annotated vOTU catalog. AMG (auxiliary
metabolic gene) prediction is a separate, more exploratory follow-on
step not included here.

## Pipeline overview
01_genomad/          geNomad virus identification on each contig file
01_phabox/           phabox2 end-to-end virus identification on each contig file
02_virus_lists/       build static per-sample fasta lists for downstream steps
03_contamination/     phabox2 contamination/provirus check (two sources)
04_checkv/             checkV quality/completeness (two sources)
05_catalog_build/      merge geNomad + phabox2 calls, trim proviruses,
dedupe, build a BLAST-safe non-redundant catalog
06_votu/               cluster into species-level vOTUs, then annotate
taxonomy / lifestyle / host prediction


Two catalogs come out of step 05:
- **main** — Complete/High/Medium CheckV quality only (the one that
  gets clustered into vOTUs in step 06)
- **relaxed** — all quality tiers included, ~100x larger, kept as a
  separate lower-confidence candidate resource

## Why geNomad AND phabox2?

Running two independently-built identification tools and taking the
union (not just the intersection) meaningfully grows the final catalog
— in practice, each tool calls a real chunk of sequences the other
misses entirely, not just overlapping calls on the same contigs. The
`called_by` column preserved through the whole pipeline (genomad /
phabox / genomad+phabox) lets you see exactly how much confidence each
final sequence has.

## Prerequisites

- SLURM cluster access
- [geNomad](https://github.com/apczok/genomad) installed, with its
  database downloaded[cite: 4]
- [phabox2](https://github.com/KennthShang/PhaBOX) installed (`conda
  activate phabox_env` or equivalent), with its database downloaded[cite: 4]
- [checkV](https://bitbucket.org/berkeleylab/checkv) installed, with
  its database downloaded[cite: 4]
- Python 3 with the packages in `requirements.txt` (`pip install -r
  requirements.txt`)[cite: 4]

## Setup

1. Copy `config.sh` and edit `BASE_DIR`, `CONTIGS_DIR`, `GENOMAD_DB`,
   `PHABOX_DB` for your environment[cite: 4].
2. In every SLURM script, replace `#SBATCH --account=CHANGE_ME` with
   your own account (and adjust `--partition` if needed)[cite: 4]. SLURM reads
   `#SBATCH` directives before any shell code runs, so these can't be
   pulled from `config.sh` — they must be edited directly in each
   script[cite: 4].
3. If your checkV module exports its database path under a different
   variable name than checkV expects (`$CHECKV_DB` vs `$CHECKVDB` is a
   common mismatch), set `export CHECKVDB=...` explicitly — see the
   commented-out line in both `04_checkv/*.sh` scripts[cite: 4].

## Running it

```bash
# 1. Identification (geNomad & PhaBOX2)
scripts/01_genomad/build_sample_manifest.sh

# Run geNomad identification
sbatch --array=1-$(wc -l < "${BASE_DIR}/contig_manifest.txt")%10 \
    scripts/01_genomad/genomad_array.sh

# Run PhaBOX2 end-to-end identification
sbatch --array=1-$(wc -l < "${BASE_DIR}/contig_manifest.txt")%10 \
    scripts/01_phabox/phabox_array.sh

# Check for incomplete samples and rerun any by name (not by index — see note below)
scripts/01_genomad/check_genomad_completeness.sh
sbatch scripts/01_genomad/genomad_single.sh /path/to/contigs/<sample>.fa

scripts/01_phabox/check_phabox_completeness.sh
sbatch scripts/01_phabox/phabox_single.sh /path/to/contigs/<sample>.fa

# 2. Build the static virus-fasta lists every later step reads from
scripts/02_virus_lists/build_virus_fasta_lists.sh

# 3. phabox2 contamination check (both sources, independent output roots)
sbatch --array=1-$(wc -l < "${BASE_DIR}/virus_fasta_list.txt")%10 \
    scripts/03_contamination/phabox_contamination_from_genomad.sh
sbatch --array=1-$(wc -l < "${BASE_DIR}/phabox_uncertain_list.txt")%10 \
    scripts/03_contamination/phabox_contamination_from_phabox.sh

# 4. checkV (both sources)
sbatch --array=1-$(wc -l < "${BASE_DIR}/virus_fasta_list.txt")%10 \
    scripts/04_checkv/checkv_from_genomad_array.sh
sbatch --array=1-$(wc -l < "${BASE_DIR}/phabox_predicted_virus_list.txt")%10 \
    scripts/04_checkv/checkv_from_phabox_array.sh

# 5. Build the non-redundant catalogs
sbatch scripts/05_catalog_build/mergetrim_main_slurm.sh
sbatch scripts/05_catalog_build/mergetrim_relaxed_slurm.sh   # optional, much larger

# 6. vOTU clustering + annotation
sbatch scripts/06_votu/phabox_votu.sh
sbatch scripts/06_votu/phabox_annotate_votu.sh
python scripts/06_votu/make_votu_table.py

Lessons baked into these scripts

A few things broke during development that are worth knowing about,
since the fixes are load-bearing, not stylistic[cite: 4]:

    Array jobs read from a static manifest file (sed -n "${SLURM_ARRAY_TASK_ID}p" manifest.txt), never from a live
    contigs/*.fa glob indexed directly by SLURM_ARRAY_TASK_ID. A
    live glob's order can shift if files are added or removed between
    computing an index and the array job actually running — this
    silently processes the wrong sample under a given index, with no
    error[cite: 4]. If you need to rerun specific failed/missing samples, rerun
    them by filename (genomad_single.sh <path> or phabox_single.sh <path>), not by resubmitting
    an array with guessed index numbers[cite: 4].

    The contamination task takes different input files depending on
    the source. For phabox2's own pipeline, the correct input is
    uncertain_sequences_for_contamination_task.fa (phabox2's own
    documented file for this task) — not predicted_virus.fa[cite: 4]. For
    geNomad's calls, it's geNomad's *_virus.fna output directly[cite: 4].

    BLAST-safe headers matter before clustering. makeblastdb -parse_seqids (used internally by phabox2 --task votu) caps local
    sequence IDs at 50 characters[cite: 4]. Long, informative headers like
    sample|contig|called_by=... will fail with "the local id is too
    long" — fix_headers.py moves that provenance into the companion TSV
    and gives the fasta short sequential IDs instead[cite: 4].

    phabox2 --task votu labels every unclustered sequence
    singleton (the literal string, not a unique ID per sequence)[cite: 4].
    make_votu_table.py gives each singleton its own ID before counting
    distinct vOTUs, or a naive .nunique() will undercount how many
    distinct viral populations you actually have[cite: 4].
