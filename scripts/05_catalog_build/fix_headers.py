#!/usr/bin/env python3
"""
Assigns short, BLAST-safe sequence IDs to a merged catalog fasta + its
companion TSV. NCBI's makeblastdb -parse_seqids caps local IDs at 50
characters; long, pipe-delimited provenance headers (sample|contig|...)
will fail there with "the local id is too long". This script keeps full
provenance in the TSV under a new `catalog_id` column instead of
encoding it in the fasta header.

Usage:
    python fix_headers.py <catalog_basename> <id_prefix>

Example:
    python fix_headers.py catalog_main_nonredundant MAIN
    -> reads    catalog_main_nonredundant.fasta / .tsv
    -> writes   catalog_main_nonredundant_safe.fasta / .tsv
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from catalog_utils import read_fasta_records  # noqa: E402


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)

    basename, prefix = sys.argv[1], sys.argv[2]
    fasta_path = Path(f"{basename}.fasta")
    tsv_path = Path(f"{basename}.tsv")

    df = pd.read_csv(tsv_path, sep="\t")
    records = read_fasta_records(fasta_path)
    assert len(records) == len(df), (
        f"mismatch: {len(records)} fasta records vs {len(df)} tsv rows — "
        "check that the fasta and tsv came from the same run"
    )

    new_ids = [f"{prefix}_{i:07d}" for i in range(1, len(df) + 1)]
    df.insert(0, "catalog_id", new_ids)

    with open(f"{basename}_safe.fasta", "w") as out:
        for new_id, (_, seq) in zip(new_ids, records):
            out.write(f">{new_id}\n{seq}\n")

    df.to_csv(f"{basename}_safe.tsv", sep="\t", index=False)
    print(f"{basename}: {len(df)} records -> {basename}_safe.fasta / .tsv")


if __name__ == "__main__":
    main()
