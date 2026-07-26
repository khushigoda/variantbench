#!/usr/bin/env python3
"""
04d — Annotation of fine-mapped variants (missense focus).

Take-home task: "Fine-mapping and analysis of missense variants shown by fine-mapping."

This step converts SuSiE's *statistical* credible sets (04c) into *biological*
mechanism. It reads only the 04c tables (never the .susie.rds, never re-fits) and
annotates every 95% credible-set variant via the Ensembl VEP REST API, then flags
the missense subset — the variants fine-mapping has nominated as causal AND that
change a protein.

Design notes
------------
* Fastest no-install VEP path: a single batch POST to /vep/human/id per chunk of
  <=200 rsIDs. All 132 CS variants are rsIDs (verified), so no allele-string build
  is needed to query. Pure stdlib (urllib + json) — runs in the `variantbench`
  (envs.main) Python with zero extra packages.
* GRCh38 endpoint (rest.ensembl.org) matches the 1000G-phase3-GRCh38 LD panel.
* Multiallelic correctness: rs11591147 (PCSK9 R46L) is G/A/T — VEP returns a
  transcript row PER alt allele (A->R/H, T->R/L). We MUST pick the row whose
  variant_allele == the GWAS effect allele (ea). The effect allele is not in the
  finemap table, so we join each CS variant to its locus gwas.tsv on (POS, z)
  — z is identical across the two files and uniquely disambiguates multiallelic
  positions where POS alone is ambiguous.
* Transcript choice: prefer MANE Select, then canonical, else most-severe.
* Raw VEP JSON is cached to <trait>.vep_raw.json so re-runs / 04e never re-hit the API.

Inputs  (env-driven, mirrors 04c):
  TRAIT     trait id (default LDL_C)
  FMDIR     dir with 04c outputs      (default results/finemap)
  LDDIR     dir with per-locus gwas   (default results/ld)
  OUTDIR    output dir                (default results/finemap)
  COVERAGE  CS coverage to annotate   (default 0.95)
  VEP_HOST  REST host                 (default https://rest.ensembl.org)
  VEP_CHUNK batch size                (default 200)

Outputs:
  <OUTDIR>/<TRAIT>.annot.tsv           one row per 95% CS variant, all consequences
  <OUTDIR>/<TRAIT>.annot.missense.tsv  the missense subset (the headline deliverable)
  <OUTDIR>/<TRAIT>.vep_raw.json        cached raw VEP response (regenerable scratch)
"""
import os, sys, csv, json, time, urllib.request, urllib.error
from collections import defaultdict

# ---------------------------------------------------------------- config
TRAIT    = os.environ.get("TRAIT",   "LDL_C")
FMDIR    = os.environ.get("FMDIR",   "results/finemap")
LDDIR    = os.environ.get("LDDIR",   "results/ld")
OUTDIR   = os.environ.get("OUTDIR",  "results/finemap")
COVERAGE = float(os.environ.get("COVERAGE", "0.95"))
VEP_HOST = os.environ.get("VEP_HOST", "https://rest.ensembl.org").rstrip("/")
VEP_CHUNK= int(os.environ.get("VEP_CHUNK", "200"))

os.makedirs(OUTDIR, exist_ok=True)
FIN   = os.path.join(FMDIR, f"{TRAIT}.finemap.tsv")
CSSUM = os.path.join(FMDIR, f"{TRAIT}.cs_summary.tsv")
RAW   = os.path.join(OUTDIR, f"{TRAIT}.vep_raw.json")
OUT   = os.path.join(OUTDIR, f"{TRAIT}.annot.tsv")
OUTM  = os.path.join(OUTDIR, f"{TRAIT}.annot.missense.tsv")

cov_tag = f"{COVERAGE:.2f}"
print(f"== 04d VEP annotation | trait={TRAIT} coverage={cov_tag} host={VEP_HOST} ==", flush=True)

# Expected effector gene(s) per locus name (locus label -> gene symbol set).
# All labels ARE the gene symbol except ABCG5_8 (a two-gene locus).
def expected_genes(locus):
    if locus == "ABCG5_8":
        return {"ABCG5", "ABCG8"}
    return {locus}

# ---------------------------------------------------------------- 1. CS variants
def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

fin = read_tsv(FIN)
cscol = "cs95" if abs(COVERAGE - 0.95) < 1e-9 else "cs99"
cs_vars = [r for r in fin if r.get(cscol, "0") not in ("0", "", "NA")]
print(f"   {len(cs_vars)} variants in {cov_tag} credible sets (col {cscol})", flush=True)

# pass_filter per (gene, cs) at this coverage, from cs_summary
passf = {}
for r in read_tsv(CSSUM):
    if abs(float(r["coverage"]) - COVERAGE) < 1e-9:
        passf[(r["gene"], r["cs"])] = (r.get("pass_filter", "NA") == "TRUE")

# ---------------------------------------------------------------- 2. effect alleles
# join each CS variant to its locus gwas.tsv on (pos, z) -> (ea, oa)
# key z to 6 dp to be robust to float formatting; pos is exact.
def zkey(z):
    return f"{float(z):.6f}"

allele = {}   # (gene, ID) -> (ea, oa)
by_locus = defaultdict(list)
for r in cs_vars:
    by_locus[r["gene"]].append(r)

for gene, rows in by_locus.items():
    gpath = os.path.join(LDDIR, f"{TRAIT}.{gene}.gwas.tsv")
    idx = {}   # (pos, zkey) -> (ea, oa)
    with open(gpath) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for g in rd:
            idx[(g["pos"], zkey(g["z"]))] = (g["ea"], g["oa"])
    for r in rows:
        allele[(gene, r["ID"])] = idx.get((r["POS"], zkey(r["z"])), ("", ""))

n_noallele = sum(1 for v in allele.values() if v == ("", ""))
print(f"   effect allele resolved for {len(allele)-n_noallele}/{len(allele)} "
      f"(missing {n_noallele})", flush=True)

# ---------------------------------------------------------------- 3. VEP (cached)
rsids = sorted({r["ID"] for r in cs_vars})
if os.path.exists(RAW):
    with open(RAW) as fh:
        vep = json.load(fh)
    print(f"   VEP: loaded cached {len(vep)} records from {os.path.basename(RAW)}", flush=True)
else:
    vep = []
    for i in range(0, len(rsids), VEP_CHUNK):
        chunk = rsids[i:i+VEP_CHUNK]
        body = json.dumps({"ids": chunk}).encode()
        req = urllib.request.Request(f"{VEP_HOST}/vep/human/id", data=body, method="POST",
                headers={"Content-Type": "application/json", "Accept": "application/json"})
        for attempt in range(4):
            try:
                with urllib.request.urlopen(req, timeout=120) as resp:
                    vep.extend(json.load(resp))
                break
            except urllib.error.HTTPError as e:
                if e.code == 429 and attempt < 3:            # rate limited: back off
                    wait = int(e.headers.get("Retry-After", 2))
                    print(f"   VEP 429, retry in {wait}s", flush=True); time.sleep(wait); continue
                raise
            except urllib.error.URLError:
                if attempt < 3:
                    time.sleep(2 ** attempt); continue
                raise
        print(f"   VEP: {min(i+VEP_CHUNK,len(rsids))}/{len(rsids)} queried", flush=True)
        time.sleep(0.34)                                     # <=3 req/s (Ensembl limit)
    with open(RAW, "w") as fh:
        json.dump(vep, fh)
    print(f"   VEP: cached {len(vep)} records -> {os.path.basename(RAW)}", flush=True)

vep_by_id = {rec.get("input", rec.get("id")): rec for rec in vep}

# ---------------------------------------------------------------- 4. pick allele-matched consequence
def pick_transcript(rec, ea):
    """From a VEP record, return the transcript_consequence for the GWAS effect
    allele, preferring MANE Select > canonical > first. Falls back across alleles
    if ea has no coding row so we still report the most severe coding change."""
    tcs = rec.get("transcript_consequences", []) or []
    # rows whose variant_allele == effect allele
    ea_rows = [t for t in tcs if t.get("variant_allele") == ea] if ea else []
    pool = ea_rows if ea_rows else tcs
    if not pool:
        return None, ("mismatch" if (ea and tcs) else "no_transcript")
    def rank(t):
        return (0 if t.get("mane_select") else 1,
                0 if t.get("canonical") else 1)
    pool = sorted(pool, key=rank)
    # prefer a coding-changing row within the pool if present
    coding = [t for t in pool if t.get("amino_acids")]
    chosen = (coding or pool)[0]
    tag = "ea" if ea_rows else ("allele_fallback" if ea else "no_ea")
    return chosen, tag

MISSENSE = "missense_variant"
KNOWN_CAUSAL = {"rs11591147"}   # PCSK9 R46L — positive control

rows_out = []
for r in cs_vars:
    gene, rid = r["gene"], r["ID"]
    ea, oa = allele.get((gene, rid), ("", ""))
    rec = vep_by_id.get(rid, {})
    msc = rec.get("most_severe_consequence", "NA")
    t, tag = pick_transcript(rec, ea)
    if t is None:
        t = {}
    cons_terms = t.get("consequence_terms", []) or []
    consequence = cons_terms[0] if cons_terms else msc
    is_missense = MISSENSE in cons_terms or (not cons_terms and msc == MISSENSE)
    aa = t.get("amino_acids", "")
    prot = t.get("protein_start", "")
    aa_change = f"{aa.split('/')[0]}{prot}{aa.split('/')[1]}" if (aa and "/" in aa and prot != "") else ""
    gsym = t.get("gene_symbol", "")
    rows_out.append({
        "trait": TRAIT, "gene": gene, "chr": r["chr"], "ID": rid, "POS": r["POS"],
        "ea": ea, "oa": oa, "pip": r["pip"], "z": r["z"], "r2_lead": r.get("r2_lead", ""),
        "cs": r[cscol], "coverage": cov_tag,
        "pass_filter": passf.get((gene, r[cscol]), ""),
        "most_severe": msc, "consequence": consequence,
        "impact": t.get("impact", ""), "is_missense": is_missense,
        "gene_symbol": gsym, "is_expected_gene": (gsym in expected_genes(gene)) if gsym else "",
        "aa_change": aa_change, "codons": t.get("codons", ""),
        "sift_pred": t.get("sift_prediction", ""), "sift_score": t.get("sift_score", ""),
        "polyphen_pred": t.get("polyphen_prediction", ""), "polyphen_score": t.get("polyphen_score", ""),
        "biotype": t.get("biotype", ""), "allele_match": tag,
        "is_known_causal": (rid in KNOWN_CAUSAL and is_missense),
    })

# sort: by chr, then pip desc within locus
rows_out.sort(key=lambda x: (int(x["chr"]), x["gene"], -float(x["pip"])))

FIELDS = ["trait","gene","chr","ID","POS","ea","oa","pip","z","r2_lead","cs","coverage",
          "pass_filter","most_severe","consequence","impact","is_missense","gene_symbol",
          "is_expected_gene","aa_change","codons","sift_pred","sift_score",
          "polyphen_pred","polyphen_score","biotype","allele_match","is_known_causal"]

def write_tsv(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)

write_tsv(OUT, rows_out)
missense = [r for r in rows_out if r["is_missense"]]
write_tsv(OUTM, missense)

# ---------------------------------------------------------------- 5. report
n_miss = len(missense)
n_miss_expected = sum(1 for r in missense if r["is_expected_gene"] is True)
loci_with_missense = sorted({r["gene"] for r in missense})
known = [r for r in rows_out if r["is_known_causal"]]
print(f"\n[done] {len(rows_out)} CS variants annotated | {n_miss} missense "
      f"({n_miss_expected} in expected effector gene) across {len(loci_with_missense)} loci",
      flush=True)
if known:
    k = known[0]
    print(f"   GROUND TRUTH: {k['ID']} {k['gene']} {k['aa_change']} "
          f"pip={k['pip']} sift={k['sift_pred']} polyphen={k['polyphen_pred']} "
          f"-> {'PASS' if k['is_expected_gene'] else 'CHECK'}", flush=True)
else:
    print("   GROUND TRUTH: rs11591147 not found as missense — INVESTIGATE", flush=True)
print(f"   -> {OUT}\n   -> {OUTM}\n   -> {RAW}", flush=True)
