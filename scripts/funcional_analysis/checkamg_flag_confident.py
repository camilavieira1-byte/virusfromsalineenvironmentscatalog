"""
Flag confident CheckAMG calls for manual review.

Reads  : <BASE>/checkamg_complete_out/03_aggregate/aggregated_results.tsv and <BASE>/complete_vOTUs.fasta
Writes : checkamg_confident_calls_flagged.tsv, checkamg_metabolic_strict.tsv,
         checkamg_denovo_only_candidates.tsv, checkamg_absent_vOTUs.tsv (in the 03_aggregate folder)

BASE defaults to the current directory; override with the AMG_BASE environment variable.
"""
import os

import pandas as pd

BASE = os.environ.get("AMG_BASE", os.getcwd())
OUT = f"{BASE}/checkamg_complete_out/03_aggregate"

# ---------- load CheckAMG aggregate table ----------
agg = pd.read_csv(f"{OUT}/aggregated_results.tsv", sep="\t")
agg["Function"] = agg["Function"].fillna("")

# ---------- vOTU lengths from the input fasta ----------
lengths, name, n = {}, None, 0
with open(f"{BASE}/complete_vOTUs.fasta") as fh:
    for line in fh:
        if line.startswith(">"):
            if name:
                lengths[name] = n
            name, n = line[1:].split()[0], 0
        else:
            n += len(line.strip())
    if name:
        lengths[name] = n

scored = set(agg["Genome"].dropna())
unscored = sorted(set(lengths) - scored)
print(f"vOTUs in fasta: {len(lengths)} | present in CheckAMG table: {len(scored)} | absent: {len(unscored)}")
if unscored:
    ul = [lengths[g] for g in unscored]
    print(f"  absent ones: {sum(l >= 5000 for l in ul)} are >=5000 bp; length range {min(ul)}-{max(ul)} bp")
    pd.DataFrame({"Genome": unscored, "length": ul}).to_csv(
        f"{OUT}/checkamg_absent_vOTUs.tsv", sep="\t", index=False)

# ---------- confident AVG calls (annotate class + high/very high) ----------
is_conf_level = agg["Final AVG Confidence Level"].isin(["high", "very high"])
conf = agg[agg["Classification"].notna() & (agg["Classification"] != "unclassified") & is_conf_level].copy()
ng = conf["Genome"].nunique()
print(f"\nconfident calls: {len(conf)} on {ng} vOTUs "
      f"({ng/len(scored)*100:.1f}% of those present, {ng/len(lengths)*100:.1f}% of all {len(lengths)})")
print(conf["Classification"].value_counts().to_string())

# ---------- review flags (tags for manual review, NOT automatic deletions) ----------
# Metabolic flags follow the caveats in Martin et al. 2025 (Nat Microbiol 10:2122)
MET_FLAGS = {
    "flag_nucleotide_core": r"ribonucleotide reductase|ribonucleoside-(?:di|tri)phosphate reductase|\bnrd[a-z]\b|dutp|dtmp kinase|thymidylate kinase|dctp deaminase",
    "flag_methyltransferase": r"(?<!hydroxy)methyltransferase|\bdcm\b|dnmt|cytosine-5",
    "flag_one_carbon": r"hydroxymethyltransferase|\bshmt\b|\bglya\b",
    "flag_glyco": r"glycosyl ?transferase|glycoside hydrolase|rhamnosyl|mannosyl|\bgt\d+\b|\bgh\d+\b",
    "flag_queuosine": r"queuosine|\bque[cdef]\b|7-cyano-7-deazaguanine|7-carboxy-7-deazaguanine|6-carboxy-5,6,7,8-tetrahydropterin|\bpreq",
    "flag_folate": r"folate|\bfol[abekp]\b|\bgch1\b|dihydrofolate|dihydroneopterin|dihydropteroate|thymidylate synthase|\bthy[ax]\b",
    "flag_adenylyltransferase": r"adenylyl ?transferase|nucleotidyl ?transferase",
}
# Regulatory flag: genes that are often part of the phage's own regulation/replication
# (heuristic keyword list; adjust to your dataset)
REG_FLAGS = {
    "flag_phage_own_regulator": r"anti-?repressor|\bcii\b|cii-like|\brha\b|rha-like|parb|spo0j|partition|gp2\.5|ssdna|single-stranded|late transcriptional activator|gapr|dnat",
}

f = conf["Function"].str.lower()
for col, pat in {**MET_FLAGS, **REG_FLAGS}.items():
    conf[col] = f.str.contains(pat, regex=True)
all_flag_cols = list(MET_FLAGS) + list(REG_FLAGS)
conf["flags"] = conf[all_flag_cols].apply(
    lambda r: ",".join(c.replace("flag_", "") for c in all_flag_cols if r[c]), axis=1)
conf.to_csv(f"{OUT}/checkamg_confident_calls_flagged.tsv", sep="\t", index=False)

# ---------- metabolic summary ----------
met = conf[conf["Classification"] == "metabolic"]
print(f"\nmetabolic confident calls: {len(met)} on {met['Genome'].nunique()} vOTUs")
for c in MET_FLAGS:
    print(f"  {c}: {int(met[c].sum())}")
strict = met[~met[list(MET_FLAGS)].any(axis=1)]
print(f"\nmetabolic, no review flag: {len(strict)} calls on {strict['Genome'].nunique()} vOTUs")
strict.to_csv(f"{OUT}/checkamg_metabolic_strict.tsv", sep="\t", index=False)
print("\ntop functions (strict metabolic):")
print(strict["Function"].value_counts().head(30).to_string())

# ---------- regulatory summary ----------
reg = conf[conf["Classification"] == "regulatory"]
print(f"\nregulatory confident calls: {len(reg)} on {reg['Genome'].nunique()} vOTUs")
print(f"  flag_phage_own_regulator: {int(reg['flag_phage_own_regulator'].sum())}")
reg_rest = reg[~reg["flag_phage_own_regulator"]]
print(f"  regulatory without that flag: {len(reg_rest)} calls on {reg_rest['Genome'].nunique()} vOTUs")
print(reg_rest["Function"].value_counts().head(15).to_string())

# ---------- de-novo-only candidates (no annotate class, high/very high) ----------
dn = agg[(agg["Classification"] == "unclassified") & is_conf_level]
dn.to_csv(f"{OUT}/checkamg_denovo_only_candidates.tsv", sep="\t", index=False)
print(f"\nde-novo-only candidates (high/very high): {len(dn)} proteins on {dn['Genome'].nunique()} vOTUs")
print(f"  of these, with an empty Function: {int((dn['Function'] == '').sum())}")
