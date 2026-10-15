"""
Relaxed genomic-context filter (Tier 2) for VIBRANT AMG candidates.

Reads  : <VIB_DIR>/vibrant_complete_tier1_FIXED.tsv and the VIBRANT annotation table
Writes : <VIB_DIR>/vibrant_complete_tier2_RELAXED.tsv

A Tier 1 gene is kept if a viral-like gene (non-zero KO/Pfam/VOG v-score) lies within 3 genes,
or if it has a viral-like direct neighbour and the genome's viral-like fraction is >= 0.20.

BASE defaults to the current directory; override with the AMG_BASE environment variable.
"""
import os
import re
import pandas as pd

BASE = os.environ.get("AMG_BASE", os.getcwd())
VIB_DIR = f"{BASE}/vibrant_amg_complete_output"
VIB_RESULTS = f"{VIB_DIR}/VIBRANT_complete_vOTUs/VIBRANT_results_complete_vOTUs"

print("[1/3] Loading VIBRANT annotation table...")
ann = pd.read_csv(f"{VIB_RESULTS}/VIBRANT_annotations_complete_vOTUs.tsv", sep="\t")

# Clean scaffold base and parse fragment/gene indices
ann["base"] = ann["scaffold"].str.replace(r"_fragment_\d+$", "", regex=True)
ann["frag"] = ann["scaffold"].str.extract(r"_fragment_(\d+)$")[0].fillna(0).astype(int)
ann["gene_no"] = ann["protein"].str.extract(r"_(\d+)$")[0].astype(int)
ann = ann.sort_values(["base", "frag", "gene_no"]).reset_index(drop=True)

# Calculate viral-like score status
vs = ann[["KO v-score", "Pfam v-score", "VOG v-score"]].apply(pd.to_numeric, errors="coerce").fillna(0)
ann["viral_like"] = (vs != 0).any(axis=1)

# Compute genomic context indicators
g = ann.groupby("base")["viral_like"]
def shift_nb(d):
    return g.shift(d, fill_value=False).astype(bool)

ann["one_neighbor"] = shift_nb(1) | shift_nb(-1)
w3 = pd.Series(False, index=ann.index)
for d in (1, 2, 3):
    w3 = w3 | shift_nb(d) | shift_nb(-d)
ann["within3"] = w3
ann["genome_viral_frac"] = ann.groupby("base")["viral_like"].transform("mean")

# Load Tier 1 dataset
t1 = pd.read_csv(f"{VIB_DIR}/vibrant_complete_tier1_FIXED.tsv", sep="\t")

# Merge context indicators into Tier 1
cols = ["protein", "base", "one_neighbor", "within3", "genome_viral_frac"]
t1 = t1.merge(ann[cols], on="protein", how="left")

# Define Relaxed Tier 2 condition: within 3 genes OR (>= 1 neighbor AND viral fraction >= 0.20)
t1["tier2_relaxed"] = t1["within3"] | (t1["one_neighbor"] & (t1["genome_viral_frac"] >= 0.20))

t2_df = t1[t1["tier2_relaxed"]].copy()

# Retain clean schema
output_cols = [c for c in t1.columns if c not in ["one_neighbor", "within3", "genome_viral_frac", "tier2_relaxed"]]
t2_df = t2_df[output_cols]

out_path = f"{VIB_DIR}/vibrant_complete_tier2_RELAXED.tsv"
t2_df.to_csv(out_path, sep="\t", index=False)

print(f"[SUMMARY] Total Tier 1 candidate genes: {len(t1)}")
print(f"[SUMMARY] Tier 2 Relaxed candidate genes retained: {len(t2_df)}")
print(f"[SUCCESS] Saved relaxed Tier 2 table to: {out_path}")
