#!/usr/bin/env python3
"""
Joins geNomad's identified viral contigs with two independent QC calls
on those same contigs — phabox2's contamination/provirus task and
checkV's quality/completeness assessment — into one per-contig table,
then adds a recommended_action column based on where the two tools
agree or disagree.

Run this AFTER:
  - 03_contamination/phabox_contamination_from_genomad.sh has completed
  - 04_checkv/checkv_from_genomad_array.sh has completed

Output: <BASE_DIR>/comparison_genomad_pipeline_annotated.tsv
        (consumed by mergetrim_main.py / mergetrim_relaxed.py)

Usage: python build_comparison_table.py
"""
import os
from pathlib import Path

import pandas as pd

BASE = Path(os.environ.get("BASE_DIR", "."))


def agreement_status(row):
    p, c = row["phabox_Provirus"], row["checkv_provirus"]
    if pd.isna(p) or pd.isna(c):
        return "not_scored_by_both"
    if p == "Yes" and c == "Yes":
        return "both_agree_provirus"
    if p == "No" and c == "No":
        return "both_agree_clean"
    if p == "Yes" and c == "No":
        return "phabox_only_provirus"
    return "checkv_only_provirus"


def host_signal_no_boundary(row):
    checkv_quality = str(row.get("checkv_checkv_quality", ""))
    checkv_host_genes = row.get("checkv_host_genes", 0) or 0
    phabox_pure_viral = str(row.get("phabox_Pure_viral", "")).lower()
    phabox_prok_genes = row.get("phabox_Prokaryotic_genes", 0) or 0

    checkv_flag = checkv_quality in ("Not-determined", "Low-quality") and checkv_host_genes > 0
    phabox_flag = "host marker" in phabox_pure_viral and phabox_prok_genes > 0
    return checkv_flag or phabox_flag


def recommended_action(row):
    status = row["agreement_status"]
    if status == "both_agree_provirus":
        return "trim_at_checkv_boundary"
    if status in ("phabox_only_provirus", "checkv_only_provirus"):
        return "flag_for_manual_review"
    if status == "both_agree_clean" and row["host_signal_no_boundary"]:
        return "flag_low_confidence_host_signal"
    return "keep_as_is"


def main():
    rows = []
    for phabox_tsv in sorted(
        (BASE / "phabox_contamination_genomad").glob("*/final_prediction/contamination_prediction.tsv")
    ):
        sample = phabox_tsv.parts[-3]
        checkv_tsv = BASE / "checkv_genomad_outputs" / sample / "quality_summary.tsv"

        pb = pd.read_csv(phabox_tsv, sep="\t").add_prefix("phabox_")
        pb = pb.rename(columns={"phabox_Accession": "contig_id"})

        if checkv_tsv.exists():
            cv = pd.read_csv(checkv_tsv, sep="\t").add_prefix("checkv_")
            cv = cv.rename(columns={"checkv_contig_id": "contig_id"})
            merged = pb.merge(cv, on="contig_id", how="outer")
        else:
            print(f"WARNING: no checkv output for {sample}")
            merged = pb

        merged.insert(0, "sample", sample)
        rows.append(merged)

    df = pd.concat(rows, ignore_index=True)

    df["agreement_status"] = df.apply(agreement_status, axis=1)
    df["host_signal_no_boundary"] = df.apply(host_signal_no_boundary, axis=1)
    df["recommended_action"] = df.apply(recommended_action, axis=1)

    out_path = BASE / "comparison_genomad_pipeline_annotated.tsv"
    df.to_csv(out_path, sep="\t", index=False)

    print(f"{df.shape[0]} contigs across {df['sample'].nunique()} samples -> {out_path}")
    print(df["recommended_action"].value_counts())


if __name__ == "__main__":
    main()
