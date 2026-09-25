#!/usr/bin/env python3
"""
Same logic as mergetrim_main.py, but with NO quality-tier filter — every
contig from both sources is included, across all five CheckV quality
tiers (Complete/High/Medium/Low-quality/Not-determined).

This is a much larger, lower-confidence catalog than the main one —
treat it as a separate candidate/supplementary resource, not something
to merge back into catalog_main_nonredundant. Most of it will be
Low-quality/Not-determined by construction.

Given the scale (expect ~100x more sequences than the main catalog),
run this as a SLURM job (see mergetrim_relaxed_slurm.sh), not
interactively.

Output: <BASE_DIR>/catalog_relaxed_nonredundant.fasta / .tsv

Usage: python mergetrim_relaxed.py
"""
import hashlib
import os
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from catalog_utils import read_fasta  # noqa: E402

BASE = Path(os.environ.get("BASE_DIR", "."))


def main():
    records = []

    # --- genomad side: no quality filter, take everything ---
    comp = pd.read_csv(BASE / "comparison_genomad_pipeline_annotated.tsv", sep="\t")

    for sample, group in comp.groupby("sample"):
        viruses = read_fasta(BASE / "checkv_genomad_outputs" / sample / "viruses.fna")
        proviruses = read_fasta(BASE / "checkv_genomad_outputs" / sample / "proviruses.fna")
        for _, row in group.iterrows():
            cid = row["contig_id"]
            if row.get("recommended_action") == "trim_at_checkv_boundary":
                seq = next((s for k, s in proviruses.items() if k.startswith(cid)), None)
            else:
                seq = viruses.get(cid)
            if seq:
                records.append(("genomad", sample, cid, seq, row["checkv_checkv_quality"]))

    # --- phabox side: no quality filter, take everything ---
    for qs_path in (BASE / "checkv_phabox_outputs").glob("*/quality_summary.tsv"):
        sample = qs_path.parts[-2]
        qs = pd.read_csv(qs_path, sep="\t")
        viruses = read_fasta(BASE / "checkv_phabox_outputs" / sample / "viruses.fna")
        proviruses = read_fasta(BASE / "checkv_phabox_outputs" / sample / "proviruses.fna")
        for _, row in qs.iterrows():
            cid = row["contig_id"]
            if row["provirus"] == "Yes":
                seq = next((s for k, s in proviruses.items() if k.startswith(cid)), None)
            else:
                seq = viruses.get(cid)
            if seq:
                records.append(("phabox", sample, cid, seq, row["checkv_quality"]))

    df = pd.DataFrame(records, columns=["origin", "sample", "contig_id", "sequence", "checkv_quality"])
    df["seq_hash"] = df["sequence"].apply(lambda s: hashlib.md5(s.encode()).hexdigest())

    before = len(df)
    origins_per_hash = df.groupby("seq_hash")["origin"].apply(lambda x: "+".join(sorted(set(x))))
    df["called_by"] = df["seq_hash"].map(origins_per_hash)
    df_dedup = df.drop_duplicates(subset="seq_hash", keep="first")
    after = len(df_dedup)

    print(f"{before} sequences before dedup, {after} after ({before - after} exact duplicates removed)")
    print(df_dedup["called_by"].value_counts())
    print(df_dedup["checkv_quality"].value_counts())

    with open(BASE / "catalog_relaxed_nonredundant.fasta", "w") as f:
        for _, row in df_dedup.iterrows():
            f.write(
                f">{row['sample']}|{row['contig_id']}|called_by={row['called_by']}"
                f"|quality={row['checkv_quality']}\n{row['sequence']}\n"
            )

    df_dedup.drop(columns="sequence").to_csv(BASE / "catalog_relaxed_nonredundant.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
