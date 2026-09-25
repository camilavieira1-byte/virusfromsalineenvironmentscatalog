#!/usr/bin/env python3
"""
Joins the per-contig catalog table (sample, origin, checkV quality —
from 05_catalog_build/fix_headers.py's *_safe.tsv output) with the vOTU
cluster assignments from phabox2 --task votu, producing one final table
that traces every original contig to its sample, origin (geNomad/
phabox/both), quality tier, and final species-level vOTU cluster.

Run this AFTER phabox_votu.sh has completed.

Output: <BASE_DIR>/votu_annotated_catalog.tsv

Usage: python make_votu_table.py
"""
import os
from pathlib import Path

import pandas as pd

BASE = Path(os.environ.get("BASE_DIR", "."))


def main():
    catalog = pd.read_csv(BASE / "catalog_main_nonredundant_safe.tsv", sep="\t")
    votu = pd.read_csv(BASE / "votu_main_catalog" / "final_prediction" / "ANI_based_vOTU.tsv", sep="\t")

    # phabox2 labels every unclustered sequence with the literal string
    # "singleton" rather than a unique ID. Left as-is, that collapses
    # every singleton into one group when you count distinct vOTUs
    # (e.g. via `.nunique()` or `groupby`) — give each singleton its own
    # ID instead, using its own accession, so cluster counts are correct.
    votu["vOTU_id"] = votu.apply(
        lambda r: r["Accession"] if r["vOTU"] == "singleton" else r["vOTU"], axis=1
    )

    merged = catalog.merge(votu, left_on="catalog_id", right_on="Accession", how="left")
    merged = merged.drop(columns=["Accession"])

    out_path = BASE / "votu_annotated_catalog.tsv"
    merged.to_csv(out_path, sep="\t", index=False)

    print(f"{merged.shape[0]} rows total -> {out_path}")
    unassigned = merged["vOTU_id"].isna().sum()
    print(f"{unassigned} contigs with no vOTU assignment (excluded by votu's default length filter)")
    print(f"{merged['vOTU_id'].nunique()} distinct vOTUs")


if __name__ == "__main__":
    main()
