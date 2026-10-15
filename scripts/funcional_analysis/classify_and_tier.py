#!/usr/bin/env python3
"""
Adds conflict classification and confidence tiers to the integrated AMG table.

Input : integrated_amg_table/integrated_amg_table_concordant.tsv   (written by build_integrated_table.py)
Output: integrated_amg_table/integrated_amg_table_tiered.tsv    (all loci + new columns)
        integrated_amg_table/integrated_amg_table_tiers1to3.tsv (tiers 1-3 only)
        integrated_amg_table/conflict_classification.tsv        (one row per conflicting tool pair)

Project directory defaults to the current directory; override with AMG_BASE (and AMG_OUT_DIR, which
must match the one used by build_integrated_table.py).

Conflict types (a conflict = two tools call the SAME physical gene but share no KO/Pfam, symbol or name):
  generic_vs_specific : one tool gives a family/domain-level name (Pfam style), the other a specific KEGG function
                        -> usually the same gene described at different levels (light review)
  related_enzyme      : both specific, same EC class (first 3 EC digits) or same enzyme type (e.g. both dehydratases)
  different_enzyme    : both specific and unrelated -> needs manual review
  generic_vs_generic  : both family-level names that do not overlap

Tiers
  T1 : all three tools support the gene, no different_enzyme conflict
  T2 : two tools support the gene and agree (shared ID/symbol/name)
  T3 : supported, but annotations disagree: two tools with a generic_vs_specific / related_enzyme / generic_vs_generic
       conflict, or three tools with a different_enzyme conflict among them, or loci linked only as ADJACENT genes
       (no coordinate overlap between tools)
  T4 : two tools, same gene, unrelated specific annotations (different_enzyme) -> review
  T5 : single tool
"""
import os
import re
import pandas as pd

BASE = os.environ.get("AMG_BASE", os.getcwd())
D = os.environ.get("AMG_OUT_DIR", f"{BASE}/integrated_amg_table")
IN = f"{D}/integrated_amg_table_concordant.tsv"

FUNC_COL = {"CheckAMG": "CheckAMG_function", "DRAM-V": "DRAMV_function", "VIBRANT": "VIBRANT_function"}
EC_RE = re.compile(r"\[EC:([^\]]+)\]")
RANK = {"different_enzyme": 3, "related_enzyme": 2, "generic_vs_generic": 1, "generic_vs_specific": 0}


def ec3(text):
    out = set()
    for blk in EC_RE.findall(text):
        for e in blk.split():
            parts = e.split(".")
            if len(parts) >= 3 and parts[0].isdigit() and parts[1].isdigit() and parts[2] != "-":
                out.add(".".join(parts[:3]))
    return out


def enzyme_types(text):
    t = EC_RE.sub(" ", text.lower())
    t = re.sub(r"\bmethylation\b", "methyltransferase", t)   # e.g. "methylation protein EvaC"
    return {w for w in re.findall(r"[a-z]+ase\b", t) if len(w) > 6}


def is_specific(desc):
    """KEGG-style description: has an EC number or starts with a gene-symbol list before ';'."""
    if EC_RE.search(desc):
        return True
    if ";" in desc:
        first = desc.split(";")[0]
        toks = [x.strip() for x in first.split(",") if x.strip()]
        return bool(toks) and all(" " not in x and len(x) <= 12 for x in toks)
    return False


def classify_pair(a, b):
    sa, sb = is_specific(a), is_specific(b)
    if sa and sb:
        if ec3(a) & ec3(b) or enzyme_types(a) & enzyme_types(b):
            return "related_enzyme"
        return "different_enzyme"
    if sa != sb:
        return "generic_vs_specific"
    return "generic_vs_generic"


def parse_basis(s):
    """'CheckAMG~DRAM-V:same_gene:ID; ...' -> [(toolA, toolB, kind, basis)]"""
    out = []
    for tok in [t.strip() for t in str(s).split(";") if t.strip()]:
        m = re.match(r"^(.+?)~(.+?):(same_gene|neighbour):(\w+)$", tok)
        if m:
            out.append(m.groups())
    return out


df = pd.read_csv(IN, sep="\t", dtype=str).fillna("")
df["n_tools"] = df["n_tools"].astype(int)

conf_rows, ctype, cpairs, tier, reason = [], [], [], [], []
for r in df.to_dict("records"):
    links = parse_basis(r["match_basis"])
    pair_types = []
    for ta, tb, kind, basis in links:
        if basis != "conflict":
            continue
        fa, fb = r.get(FUNC_COL[ta], ""), r.get(FUNC_COL[tb], "")
        t = classify_pair(fa, fb)
        pair_types.append((ta, tb, t))
        conf_rows.append({"group_id": r["group_id"], "vOTU": r["vOTU"], "tool_a": ta, "function_a": fa,
                          "tool_b": tb, "function_b": fb, "conflict_type": t})
    worst = max((t for _, _, t in pair_types), key=RANK.get, default="")
    ctype.append(worst)
    cpairs.append("; ".join(f"{a}~{b}:{t}" for a, b, t in pair_types))

    n = r["n_tools"]
    agree = [bs for _, _, _, bs in links if bs != "conflict"]
    has_same_gene = any(k == "same_gene" for _, _, k, _ in links)
    if n >= 2 and not has_same_gene:
        tier.append("T3"); reason.append("adjacent genes with the same function, not the same gene")
    elif n == 3:
        if worst == "different_enzyme":
            tier.append("T3"); reason.append("3 tools, but an unrelated-enzyme conflict among them")
        else:
            tier.append("T1"); reason.append("3 tools" + (f", conflict only {worst}" if worst else ""))
    elif n == 2:
        if worst == "":
            how = "ID/symbol" if any(x in ("ID", "symbol") for x in agree) else "name"
            tier.append("T2"); reason.append(f"2 tools agree ({how})")
        elif worst == "different_enzyme":
            tier.append("T4"); reason.append("2 tools, same gene, unrelated specific annotations")
        else:
            tier.append("T3"); reason.append(f"2 tools, same gene, {worst}")
    else:
        tier.append("T5"); reason.append("single tool")

df["conflict_type"] = ctype
df["conflict_pairs"] = cpairs
df["confidence_tier"] = tier
df["tier_reason"] = reason

out_all = f"{D}/integrated_amg_table_tiered.tsv"
df.to_csv(out_all, sep="\t", index=False)
df[df["confidence_tier"].isin(["T1", "T2", "T3"])].to_csv(f"{D}/integrated_amg_table_tiers1to3.tsv", sep="\t", index=False)
pd.DataFrame(conf_rows).to_csv(f"{D}/conflict_classification.tsv", sep="\t", index=False)

print(f"loci: {len(df)}")
print(df["confidence_tier"].value_counts().sort_index().to_string())
print("\ntier x tools_supporting:")
print(pd.crosstab(df["confidence_tier"], df["tools_supporting"]).to_string())
print("\nconflict pairs by type:")
print(pd.DataFrame(conf_rows)["conflict_type"].value_counts().to_string() if conf_rows else "(none)")
print("\nloci per tier on vOTUs:", df.groupby("confidence_tier")["vOTU"].nunique().to_dict())
print(f"\nsaved: {out_all}")
if conf_rows:
    print("\nfirst different_enzyme / related_enzyme pairs to review:")
    c = pd.DataFrame(conf_rows)
    print(c[c.conflict_type.isin(["different_enzyme", "related_enzyme"])][
        ["group_id", "tool_a", "function_a", "tool_b", "function_b", "conflict_type"]].head(30).to_string(index=False, max_colwidth=45))
