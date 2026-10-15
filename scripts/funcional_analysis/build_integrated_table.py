#!/usr/bin/env python3
"""
Integrated AMG table: CheckAMG (annotate) + DRAM-v tier 2 + VIBRANT tier 2.

Genes from different tools are merged into one locus only if they are close in gene
index AND agree in function (shared KO/Pfam, shared gene symbol, or matching name).
Positional neighbours that disagree are NOT merged and are written to a QC file.
"""
import importlib.util
import os
import re
from collections import defaultdict
import pandas as pd

# ====================== CONFIG ======================
BASE = os.environ.get("AMG_BASE", os.getcwd())   # project directory (default: current directory)
CK = f"{BASE}/checkamg_complete_out"
AGG = f"{CK}/03_aggregate"
RES = f"{CK}/01_annotate/results"


def _checkamg_files_dir():
    """Locate CheckAMG's bundled reference tables (AMGs.tsv, AMGs_categorized.tsv)."""
    spec = importlib.util.find_spec("CheckAMG")
    if spec and spec.submodule_search_locations:
        return os.path.join(list(spec.submodule_search_locations)[0], "files")
    return ""


PKG = os.environ.get("CHECKAMG_FILES", _checkamg_files_dir())
DRAMV_FILE = f"{BASE}/dramv_distill_results/dramv_complete_tier2_FIXED.tsv"
VIBRANT_FILE = f"{BASE}/vibrant_amg_complete_output/vibrant_complete_tier2_FRAGMENTFIXED.tsv"
VIB_PATHWAYS = (f"{BASE}/vibrant_amg_complete_output/VIBRANT_complete_vOTUs/"
                f"VIBRANT_results_complete_vOTUs/VIBRANT_AMG_pathways_complete_vOTUs.tsv")
FASTA = f"{BASE}/complete_vOTUs.fasta"
OUT_DIR = os.environ.get("AMG_OUT_DIR", f"{BASE}/integrated_amg_table")
# Regex that extracts the vOTU id from scaffold/protein names (dataset-specific)
VOTU_ID_REGEX = os.environ.get("VOTU_ID_REGEX", r"(MAIN_\d+)")

VIRAL_CONF_KEEP = ("high", "medium")   # CheckAMG annotate viral-origin confidence
IDX_TOL_NAME = 1                       # max gene-index offset when matched by name only
IDX_TOL_ID = 2                         # max gene-index offset when matched by KO/Pfam or gene symbol
# Genes carrying any of these review flags are REMOVED before integration (set() = remove nothing).
# Available flags: ribosomal, peptidase, glycosyltransferase_hydrolase, methyltransferase,
#                  adenylyl_nucleotidyltransferase, folate, queuosine, nucleotide_core
EXCLUDE_FLAGS = {"ribosomal", "peptidase"}
# =====================================================

os.makedirs(OUT_DIR, exist_ok=True)

KO_RE = re.compile(r"\bK\d{5}\b")
PF_RE = re.compile(r"\bPF\d{5}\b")
CORE_RE = (r"ribonucleotide reductase|ribonucleoside-(?:di|tri)phosphate reductase|\bnrd[a-z]\b|dutp|"
           r"dtmp kinase|thymidylate kinase|dctp deaminase|queuosine|\bque[cdef]\b")
DB_COLS = ["KEGG_hmm_id", "FOAM_hmm_id", "Pfam_hmm_id", "dbCAN_hmm_id", "METABOLIC_hmm_id", "CAMPER_hmm_id"]
TOOLS = ("CheckAMG", "DRAM-V", "VIBRANT")

FLAG_RE = {k: re.compile(v) for k, v in {
    "ribosomal": r"ribosom|\brp[sl][a-z0-9]{1,2}\b",
    "peptidase": r"peptidase|protease|proteinase|\blapg\b",
    "glycosyltransferase_hydrolase": r"glycosyl ?transferase|glycoside hydrolase|rhamnosyl|mannosyl|\bgt\d+\b|\bgh\d+\b",
    "methyltransferase": r"(?<!hydroxy)methyltransferase|\bdcm\b|\bdnmt|cytosine-5",
    "adenylyl_nucleotidyltransferase": r"adenylyl ?transferase|nucleotidyl ?transferase|cytidylyl ?transferase",
    "folate": r"folate|\bfol[abekp]\b|\bgch1\b|dihydrofolate|dihydroneopterin|dihydropteroate|thymidylate synthase",
    "queuosine": r"queuosine|\bque[cdef]\b|deazaguanine",
    "nucleotide_core": (r"ribonucleotide reductase|ribonucleoside-(?:di|tri)phosphate reductase|\bnrd[a-z]\b|"
                        r"dutp|dtmp kinase|thymidylate kinase|dctp deaminase"),
}.items()}

GENERIC = set("""family domain protein proteins like type subunit putative probable from enzyme terminal class
superfamily binding catalytic containing associated related dependent component the and with for complex chain alpha
beta gamma small large short long dehydrogenase reductase oxidase oxygenase hydrolase kinase synthase synthetase
transferase isomerase epimerase dehydratase ligase lyase hydroxylase phosphatase transporter regulator
oxidoreductase motif repeat central middle region fold domains uncharacterized unknown hypothetical""".split())


# ====================== helpers ======================
def rd(path, **kw):
    return pd.read_csv(path, sep="\t", dtype=str, low_memory=False, **kw).fillna("")


def need(df, cols, name):
    miss = [c for c in cols if c not in df.columns]
    if miss:
        raise SystemExit(f"{name}: missing columns {miss}\navailable: {list(df.columns)}")


def ids_in(*texts):
    out = set()
    for t in texts:
        t = "" if t is None else str(t)
        out |= {"KO:" + k for k in KO_RE.findall(t)} | {"PF:" + p for p in PF_RE.findall(t)}
    return out


def genome_of(x):
    m = re.search(VOTU_ID_REGEX, str(x))
    return m.group(1) if m else None


def idx_of(x):
    m = re.search(r"_(\d+)$", str(x))
    return int(m.group(1)) if m else None


def join_u(vals, sep=" | "):
    seen, out = set(), []
    for v in vals:
        v = str(v).strip()
        if v and v not in seen:
            seen.add(v)
            out.append(v)
    return sep.join(out)


def fmax(vals):
    xs = []
    for v in vals:
        try:
            xs.append(float(v))
        except (TypeError, ValueError):
            pass
    return round(max(xs), 4) if xs else ""


def flags_of(text):
    t = str(text).lower()
    return {k for k, r in FLAG_RE.items() if r.search(t)}


def split_desc(desc):
    """KEGG-style 'sym, SYM; description [EC:..]' -> (gene symbols, description text)."""
    s = re.sub(r"\[[^\]]*\]", " ", str(desc).lower())
    segs = [x.strip() for x in s.split(";") if x.strip()]
    symbols, names = set(), []
    for k, seg in enumerate(segs):
        parts = [p.strip() for p in seg.split(",") if p.strip()]
        # a symbol list is never the last segment (the last one is always the description)
        is_sym = (len(segs) > 1 and k < len(segs) - 1 and parts
                  and all(" " not in p and len(p) <= 12 for p in parts))
        if is_sym:
            symbols |= {p for p in parts if len(p) >= 3 and not re.match(r"^(e?\d|ec)", p)}
        else:
            names.append(seg)
    return symbols, " ".join(names)


def name_terms(text):
    s = re.sub(r"(?<=[a-z])-(?=[a-z])", "", text.lower())
    return {t for t in re.findall(r"[a-z0-9]+", s) if len(t) >= 3 and not t.isdigit() and t not in GENERIC}


def make_node(tool, genome, gene, ids, descs, extra):
    sym, terms = set(), set()
    for dsc in descs:
        s, nm = split_desc(dsc)
        sym |= s
        terms |= name_terms(nm)
    text = " ".join(descs) + " " + extra.get("flag_text", "")
    node = {"tool": tool, "genome": genome, "gene": gene, "idx": idx_of(gene), "ids": set(ids),
            "descs": [x for x in descs if x], "symbols": sym, "terms": terms, "flags": flags_of(text)}
    node.update(extra)
    return node


def concordant(a, b):
    if a["ids"] & b["ids"]:
        return "ID"
    if a["symbols"] & b["symbols"]:
        return "symbol"
    ta, tb = a["terms"], b["terms"]
    if ta and tb:
        inter = ta & tb
        if inter and (len(inter) / len(ta | tb) >= 0.5 or ta <= tb or tb <= ta):
            return "name"
    return None


# ====================== vOTU universe ======================
lengths, name, ln = {}, None, 0
with open(FASTA) as fh:
    for line in fh:
        if line.startswith(">"):
            if name:
                lengths[name] = ln
            name, ln = line[1:].split()[0], 0
        else:
            ln += len(line.strip())
    if name:
        lengths[name] = ln
scored = set(rd(f"{AGG}/aggregated_results.tsv", usecols=["Genome"])["Genome"])
print(f"vOTUs in fasta: {len(lengths)} | analysed by CheckAMG: {len(scored)}")

# ====================== CheckAMG reference tables ======================
W, LVL, CAT = {}, {}, {}
p = f"{PKG}/AMGs.tsv"
if os.path.exists(p):
    a = rd(p)
    need(a, ["id", "AMG_weight", "amg_level"], "AMGs.tsv")
    for i, wt, lv in zip(a["id"], pd.to_numeric(a["AMG_weight"], errors="coerce"), a["amg_level"]):
        if pd.isna(wt):
            continue
        for key in {i, i.split(".")[0]}:
            if key not in W or wt > W[key]:
                W[key], LVL[key] = float(wt), lv
    print(f"AMGs.tsv: {len(W)} lookup keys")
else:
    print("WARNING: AMGs.tsv not found; AMG weights will be empty:", p)

p = f"{PKG}/AMGs_categorized.tsv"
if os.path.exists(p):
    c = rd(p)
    idc = "id" if "id" in c.columns else c.columns[0]
    cat_cols = [x for x in c.columns if x.lower().startswith("category")]
    print("AMGs_categorized.tsv columns:", list(c.columns), "-> using", cat_cols[:1])
    if cat_cols:
        tmp = defaultdict(set)
        for i, v in zip(c[idc], c[cat_cols[0]]):
            if v.strip():
                tmp[i].add(v.strip())
                tmp[i.split(".")[0]].add(v.strip())
        CAT = {k: "; ".join(sorted(v)) for k, v in tmp.items()}
else:
    print("WARNING: AMGs_categorized.tsv not found; CheckAMG categories will be empty")

VIB_CAT = {}
if os.path.exists(VIB_PATHWAYS):
    vp = rd(VIB_PATHWAYS)
    if "Present AMG KOs" in vp.columns and "Metabolism" in vp.columns:
        tmp = defaultdict(set)
        for kos, met_ in zip(vp["Present AMG KOs"], vp["Metabolism"]):
            for k in KO_RE.findall(kos):
                tmp[k].add(met_.strip())
        VIB_CAT = {k: "; ".join(sorted(v)) for k, v in tmp.items()}
        print(f"VIBRANT pathway map: {len(VIB_CAT)} KOs")


def best_weight(row):
    """Highest AMG weight over ALL database hits of a protein (KEGG, FOAM, Pfam, dbCAN, METABOLIC, CAMPER)."""
    best = ("", "", -1.0)
    for col in DB_COLS:
        i = row[col]
        if not i:
            continue
        for key in (i, i.split(".")[0]):
            if key in W and W[key] > best[2]:
                best = (i, col[:-7], W[key])
    return best


# ====================== CheckAMG nodes ======================
d = rd(f"{AGG}/aggregated_results_detailed.tsv")
need(d, ["Protein", "Genome", "Classification (annotate)", "Function (annotate)",
         "Viral Confidence Level (annotate)", "Final AVG Probability (de-novo)",
         "Final AVG Confidence Level (de-novo)"], "aggregated_results_detailed.tsv")
fr = rd(f"{RES}/final_results.tsv", usecols=["Protein", "KEGG KO", "Pfam Accession"]).drop_duplicates("Protein")
ga = rd(f"{RES}/gene_annotations.tsv", usecols=["Protein"] + DB_COLS).drop_duplicates("Protein")
met = d[d["Classification (annotate)"] == "metabolic"].merge(fr, on="Protein", how="left").merge(ga, on="Protein", how="left").fillna("")
met["ids"] = [ids_in(k, pf, kh, ph) for k, pf, kh, ph in
              zip(met["KEGG KO"], met["Pfam Accession"], met["KEGG_hmm_id"], met["Pfam_hmm_id"])]
conf_ok = met["Viral Confidence Level (annotate)"].isin(VIRAL_CONF_KEEP)
core = met["Function (annotate)"].str.lower().str.contains(CORE_RE) | met["ids"].apply(lambda s: "KO:K00558" in s)
sel = met[conf_ok & ~core]
print(f"CheckAMG metabolic proteins: {len(met)} -> viral confidence {VIRAL_CONF_KEEP}: {int(conf_ok.sum())} "
      f"-> after core-gene exclusion: {len(sel)} proteins on {sel['Genome'].nunique()} vOTUs")

raw_nodes = []
for r in sel.to_dict("records"):
    wid, wdb, wval = best_weight(r)
    cat = CAT.get(wid, CAT.get(wid.split(".")[0], "")) if wid else ""
    raw_nodes.append(make_node("CheckAMG", r["Genome"], r["Protein"], r["ids"], [r["Function (annotate)"]], {
        "viral_conf": r["Viral Confidence Level (annotate)"],
        "dn_prob": r["Final AVG Probability (de-novo)"], "dn_conf": r["Final AVG Confidence Level (de-novo)"],
        "consensus": r["Final AVG Confidence Level (de-novo)"] in ("high", "very high"),
        "w": round(wval, 3) if wval >= 0 else "", "w_id": wid, "w_db": wdb,
        "w_level": LVL.get(wid, LVL.get(wid.split(".")[0], "")) if wid else "", "cat": cat}))

# ====================== DRAM-v nodes ======================
dv_rows = rd(DRAMV_FILE)
need(dv_rows, ["gene", "gene_id", "gene_description"], "DRAM-v tier 2")
acc = {}
for r in dv_rows.to_dict("records"):
    g = genome_of(r.get("orig_scaffold") or r.get("scaffold"))
    if g is None:
        continue
    e = acc.setdefault((g, r["gene"]), {"ids": set(), "descs": [], "aux": [], "dflags": [], "hdr": [], "sub": [], "mod": []})
    e["ids"] |= ids_in(r["gene_id"])
    e["descs"].append(r["gene_description"])
    e["aux"].append(r.get("auxiliary_score", ""))
    e["dflags"].append(r.get("amg_flags", ""))
    e["hdr"].append(r.get("header", ""))
    e["sub"].append(r.get("subheader", ""))
    e["mod"].append(r.get("module", ""))
for (g, gene), e in acc.items():
    raw_nodes.append(make_node("DRAM-V", g, gene, e["ids"], e["descs"], {
        "aux": fmax(e["aux"]), "dflags": join_u(e["dflags"], ","), "hdr": join_u(e["hdr"]), "sub": join_u(e["sub"]),
        "flag_text": " ".join(e["hdr"] + e["sub"] + e["mod"])}))

# ====================== VIBRANT nodes ======================
vb_rows = rd(VIBRANT_FILE)
need(vb_rows, ["protein", "scaffold", "AMG KO", "AMG KO name", "Pfam", "Pfam name"], "VIBRANT tier 2")
acc = {}
for r in vb_rows.to_dict("records"):
    g = genome_of(r["scaffold"])
    if g is None:
        continue
    e = acc.setdefault((g, r["protein"]), {"ids": set(), "descs": [], "kos": set()})
    e["ids"] |= ids_in(r["AMG KO"], r["Pfam"])
    e["kos"] |= set(KO_RE.findall(r["AMG KO"]))
    e["descs"].append(r["AMG KO name"] or r["Pfam name"])
for (g, gene), e in acc.items():
    cats = sorted({VIB_CAT[k] for k in e["kos"] if k in VIB_CAT})
    raw_nodes.append(make_node("VIBRANT", g, gene, e["ids"], e["descs"], {"cat": "; ".join(cats)}))

# ====================== review-flag exclusion ======================
nodes, excluded = [], []
for n in raw_nodes:
    if n["flags"] & EXCLUDE_FLAGS:
        excluded.append({"tool": n["tool"], "vOTU": n["genome"], "gene": n["gene"],
                         "function": join_u(n["descs"]), "flags": ",".join(sorted(n["flags"]))})
    else:
        nodes.append(n)
print("\ngenes per tool before exclusion:", {t: sum(n["tool"] == t for n in raw_nodes) for t in TOOLS})
print(f"removed by EXCLUDE_FLAGS {sorted(EXCLUDE_FLAGS)}:", {t: sum(x['tool'] == t for x in excluded) for t in TOOLS})
print("genes per tool after exclusion: ", {t: sum(n["tool"] == t for n in nodes) for t in TOOLS})
pd.DataFrame(excluded).to_csv(f"{OUT_DIR}/excluded_genes.tsv", sep="\t", index=False)

# ====================== link genes across tools ======================
parent = list(range(len(nodes)))


def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb:
        parent[rb] = ra


by_votu = defaultdict(list)
for u, n in enumerate(nodes):
    if n["idx"] is not None:
        by_votu[n["genome"]].append(u)

link_basis, rejected = {}, []
for votu, us in by_votu.items():
    us.sort(key=lambda u: nodes[u]["idx"])
    for pos, u in enumerate(us):
        for v in us[pos + 1:]:
            na, nb = nodes[u], nodes[v]
            di = nb["idx"] - na["idx"]
            if di > max(IDX_TOL_ID, IDX_TOL_NAME):
                break
            if na["tool"] == nb["tool"]:
                continue
            basis = concordant(na, nb)
            if (basis in ("ID", "symbol") and di <= IDX_TOL_ID) or (basis == "name" and di <= IDX_TOL_NAME):
                union(u, v)
                link_basis[(u, v)] = basis
            elif di <= IDX_TOL_NAME:
                rejected.append({"vOTU": votu, "tool_a": na["tool"], "gene_a": na["gene"], "function_a": join_u(na["descs"]),
                                 "tool_b": nb["tool"], "gene_b": nb["gene"], "function_b": join_u(nb["descs"])})
pd.DataFrame(rejected).to_csv(f"{OUT_DIR}/rejected_positional_pairs.tsv", sep="\t", index=False)

# ====================== one row per locus ======================
groups = defaultdict(list)
for u in range(len(nodes)):
    groups[find(u)].append(u)
order = {t: i for i, t in enumerate(TOOLS)}
bases_by_root = defaultdict(set)
for (u, v), b in link_basis.items():
    ta, tb = sorted([nodes[u]["tool"], nodes[v]["tool"]], key=order.get)
    bases_by_root[find(u)].add(f"{ta}~{tb}:{b}")

rows = []
for root, us in groups.items():
    members = [nodes[u] for u in us]
    by = {t: [n for n in members if n["tool"] == t] for t in TOOLS}
    c, dv, vb = by["CheckAMG"], by["DRAM-V"], by["VIBRANT"]
    genome = members[0]["genome"]
    tools = [t for t in TOOLS if by[t]]
    allids = set().union(*[n["ids"] for n in members])
    if c:
        func, src = join_u(x for n in c for x in n["descs"]), "CheckAMG"
    elif dv:
        func, src = join_u(x for n in dv for x in n["descs"]), "DRAM-V"
    else:
        func, src = join_u(x for n in vb for x in n["descs"]), "VIBRANT"
    idxs = lambda ns: ",".join(map(str, sorted({n["idx"] for n in ns if n["idx"] is not None})))
    rows.append({
        "vOTU": genome, "vOTU_length": lengths.get(genome, ""), "vOTU_in_CheckAMG_run": genome in scored,
        "tools_supporting": "+".join(tools), "n_tools": len(tools),
        "function": func, "function_source": src,
        "KOs": ",".join(sorted(i[3:] for i in allids if i.startswith("KO:"))),
        "Pfams": ",".join(sorted(i[3:] for i in allids if i.startswith("PF:"))),
        "match_basis": "; ".join(sorted(bases_by_root[root])),
        "review_flags": ",".join(sorted(set().union(*[n["flags"] for n in members]))),
        "CheckAMG_proteins": join_u((n["gene"] for n in c), ","), "idx_CheckAMG": idxs(c),
        "CheckAMG_function": join_u(x for n in c for x in n["descs"]),
        "CheckAMG_viral_conf_annotate": join_u(n["viral_conf"] for n in c),
        "CheckAMG_consensus_denovo": any(n["consensus"] for n in c),
        "CheckAMG_denovo_AVG_prob": fmax(n["dn_prob"] for n in c),
        "CheckAMG_AMG_weight": fmax(n["w"] for n in c), "CheckAMG_weight_db": join_u((n["w_db"] for n in c), ","),
        "CheckAMG_AMG_level": join_u(n["w_level"] for n in c), "CheckAMG_category": join_u(n["cat"] for n in c),
        "DRAMV_genes": join_u((n["gene"] for n in dv), ","), "idx_DRAMV": idxs(dv),
        "DRAMV_function": join_u(x for n in dv for x in n["descs"]),
        "DRAMV_aux_score": fmax(n["aux"] for n in dv), "DRAMV_flags": join_u((n["dflags"] for n in dv), ","),
        "DRAMV_header": join_u(n["hdr"] for n in dv), "DRAMV_subheader": join_u(n["sub"] for n in dv),
        "VIBRANT_proteins": join_u((n["gene"] for n in vb), ","), "idx_VIBRANT": idxs(vb),
        "VIBRANT_function": join_u(x for n in vb for x in n["descs"]), "VIBRANT_category": join_u(n["cat"] for n in vb),
        "_first_idx": min((n["idx"] for n in members if n["idx"] is not None), default=0),
    })

df = pd.DataFrame(rows).sort_values(["n_tools", "vOTU", "_first_idx"], ascending=[False, True, True]).reset_index(drop=True)
df = df.drop(columns="_first_idx")
df.insert(0, "group_id", [f"L{i + 1:04d}" for i in range(len(df))])
df.to_csv(f"{OUT_DIR}/integrated_amg_table_concordant.tsv", sep="\t", index=False)
df[df["n_tools"] >= 2].to_csv(f"{OUT_DIR}/integrated_amg_table_concordant_2plus_tools.tsv", sep="\t", index=False)

# ====================== summary ======================
print(f"\nloci: {len(df)} on {df['vOTU'].nunique()} vOTUs")
print(df["tools_supporting"].value_counts().to_string())
for t in TOOLS:
    s = df[df["tools_supporting"].str.contains(t, regex=False)]
    print(f"  {t}: {len(s)} loci on {s['vOTU'].nunique()} vOTUs")
print(f"vOTUs with a locus supported by >=2 tools: {df[df.n_tools >= 2]['vOTU'].nunique()} | by all 3: {df[df.n_tools == 3]['vOTU'].nunique()}")
print("\nlinks kept, by tool pair and basis:")
cnt = defaultdict(int)
for (u, v), b in link_basis.items():
    ta, tb = sorted([nodes[u]["tool"], nodes[v]["tool"]], key=order.get)
    cnt[f"{ta}~{tb}:{b}"] += 1
for k in sorted(cnt):
    print(f"  {k}: {cnt[k]}")
print(f"positional neighbours rejected for disagreeing in function: {len(rejected)} (see rejected_positional_pairs.tsv)")
ck = df[df["CheckAMG_proteins"] != ""]
wv = pd.to_numeric(ck["CheckAMG_AMG_weight"], errors="coerce")
print(f"\nCheckAMG loci: {len(ck)} | with weight: {int(wv.notna().sum())} | weight < 0.6: {int((wv < 0.6).sum())}")
print("loci with review flags:", int((df["review_flags"] != "").sum()))
print(df["review_flags"].replace("", "(none)").str.split(",").explode().value_counts().head(10).to_string())
print("\nloci supported by all three tools:")
print(df[df.n_tools == 3][["group_id", "vOTU", "function", "match_basis"]].head(40).to_string(index=False, max_colwidth=50))
print(f"\nsaved: {OUT_DIR}/integrated_amg_table_concordant.tsv (+ _2plus_tools, rejected_positional_pairs, excluded_genes)")
