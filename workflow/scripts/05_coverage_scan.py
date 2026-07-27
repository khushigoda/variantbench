#!/usr/bin/env python3
# =============================================================================
# 05_coverage_scan.py  —  which eQTL Catalogue datasets fine-map our target genes?
# -----------------------------------------------------------------------------
# WHY THIS EXISTS
#   Liver is the biologically-obvious tissue for LDL, but GTEx liver (QTD000266,
#   n=208) is the ONLY liver expression dataset in the eQTL Catalogue and at that
#   sample size it fine-maps none of our 12 LDL effector genes. The authors (and
#   Open Targets / GTEx coloc papers) never rely on one tissue: they scan the
#   WHOLE catalogue and take the best coloc per gene. This script reproduces that
#   logic empirically — across EVERY dataset and EVERY quantification method in
#   the Catalogue (ge, exon, tx, txrev, leafcutter, microarray, aptamer), it asks
#   which datasets actually have a fine-mapped credible set for each effector
#   gene, and how well-powered (sample size, max PIP) each is.
#
#   Scanning ALL quant methods — not ge alone — is deliberate: a gene like HMGCR
#   acts on LDL through a splice-altering variant that is invisible to a total-
#   expression (ge) scan but visible to leafcutter / txrev (splicing) QTLs. A
#   ge-only coverage scan would silently drop those loci.
#
#   It reads only the small per-dataset credible_sets.tsv.gz files (~1 MB each),
#   NOT the 500 MB lbf_variable matrices — so it is light and safe to run.
#   Requests are paced (default 0.2 s) so the EBI FTP firewall does not flag us.
#
# OUTPUT (both under results/coloc/):
#   <TRAIT>.eqtl_coverage.tsv   one row per (gene x dataset) that has a CS:
#                               gene, ensg, dataset, study, tissue, quant, n, n_cs, max_pip
#   <TRAIT>.eqtl_coverage.best.tsv  one row per (gene x quant_method) = the
#                               highest-N dataset that fine-maps that gene under
#                               that quant method. Keeping the split per quant
#                               method means a splice-QTL hit (leafcutter) is not
#                               masked by a better-powered ge dataset with no
#                               shared causal variant. This is what 05a colocs.
#
# USAGE
#   TRAIT=LDL_C \
#   TARGETS=config/coloc_targets.LDL_C.tsv \
#   OUTDIR=results/coloc \
#   python workflow/scripts/05_coverage_scan.py
#
#   Optional env:  QUANT=all   SLEEP=0.2   API=https://www.ebi.ac.uk/eqtl/api/v2
#                  FTP=https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie
#                  QUANT restricts the sweep to one method (ge/exon/tx/txrev/
#                  leafcutter/microarray/aptamer); default "all" = whole catalogue.
#                  MAXDS=N caps datasets scanned (smoke tests only).
# =============================================================================
import os, sys, json, gzip, time, urllib.request, urllib.error

TRAIT   = os.environ.get("TRAIT", "LDL_C")
TARGETS = os.environ.get("TARGETS", f"config/coloc_targets.{TRAIT}.tsv")
OUTDIR  = os.environ.get("OUTDIR", "results/coloc")
QUANT   = os.environ.get("QUANT", "all")         # "all" = every quant method
SLEEP   = float(os.environ.get("SLEEP", "0.2"))  # pace FTP requests (DoS-safe)
MAXDS   = int(os.environ.get("MAXDS", "0"))      # 0 = no cap (smoke tests only)
API     = os.environ.get("API", "https://www.ebi.ac.uk/eqtl/api/v2")
FTP     = os.environ.get("FTP", "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie")

os.makedirs(OUTDIR, exist_ok=True)

def get_json(url, timeout=60):
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return json.load(r)

# ---- 1. load target genes (locus<TAB>ensg<TAB>symbol) ----------------------
targets = {}   # ensg -> symbol
with open(TARGETS) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    ie, isym = hdr.index("ensg"), hdr.index("symbol")
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) > max(ie, isym):
            targets[f[ie]] = f[isym]
print(f"== 05 coverage scan | trait={TRAIT} | {len(targets)} target genes | quant={QUANT} ==",
      flush=True)

# ---- 2. enumerate datasets (all quant methods unless QUANT restricts) ------
qfilter = "" if QUANT == "all" else f"quant_method={QUANT}&"
datasets, start = [], 0
while True:
    page = get_json(f"{API}/datasets?{qfilter}size=100&start={start}")
    if not page:
        break
    datasets += page
    if len(page) < 100:
        break
    start += 100
if MAXDS:
    datasets = datasets[:MAXDS]
from collections import Counter
qc = Counter(x.get("quant_method") for x in datasets)
print(f"   {len(datasets)} dataset(s) to scan  ({dict(qc)})\n", flush=True)

# ---- 3. scan each dataset's small credible_sets file -----------------------
# coverage row keyed by (sym, quant): we keep the per-quant-method split so a
# splice hit is not collapsed into a better-powered ge dataset downstream.
cover = []          # list of dataset-level dicts (one per gene x dataset w/ a CS)
scanned = missing = 0
for i, x in enumerate(datasets, 1):
    qtd, qts = x["dataset_id"], x["study_id"]
    tissue, study, n = x["tissue_label"], x["study_label"], x["sample_size"]
    quant = x.get("quant_method", "NA")
    url = f"{FTP}/{qts}/{qtd}/{qtd}.credible_sets.tsv.gz"
    try:
        raw = urllib.request.urlopen(url, timeout=90).read()
        txt = gzip.decompress(raw).decode()
        scanned += 1
    except Exception:
        missing += 1
        continue
    lines = txt.split("\n")
    hdr = lines[0].split("\t")
    # gene_id maps the molecular trait (exon/tx/splice ids) back to the gene.
    gi = hdr.index("gene_id") if "gene_id" in hdr else hdr.index("molecular_trait_id")
    mi = hdr.index("molecular_trait_id")
    pi = hdr.index("pip")
    agg = {}   # ensg -> [n_cs_rows, max_pip, set(mtids), best_mtid_at_max_pip]
    for line in lines[1:]:
        if not line:
            continue
        f = line.split("\t")
        g = f[gi]
        if g in targets:
            try:
                p = float(f[pi])
            except (ValueError, IndexError):
                p = 0.0
            d = agg.setdefault(g, [0, 0.0, set(), "NA"])
            d[0] += 1
            if p > d[1]:               # remember which molecular trait carries the top PIP
                d[1] = p
                d[3] = f[mi]
            d[2].add(f[mi])
    for ensg, (ncs, mp, mtids, best_mtid) in agg.items():
        cover.append(dict(gene=targets[ensg], ensg=ensg, dataset=qtd, study_id=qts,
                          study=study, tissue=tissue, quant=quant,
                          n=int(n), n_cs=ncs, max_pip=round(mp, 4),
                          n_mtid=len(mtids), best_mtid=best_mtid))
    if i % 25 == 0:
        print(f"   ...{i}/{len(datasets)} scanned  (hits so far: {len(cover)})", flush=True)
    time.sleep(SLEEP)

print(f"\n   scanned={scanned}  unreachable={missing}  coverage rows={len(cover)}", flush=True)

# ---- 4. write full coverage + best-per-(gene x quant) tables ---------------
cols = ["gene", "ensg", "dataset", "study_id", "study", "tissue", "quant",
        "n", "n_cs", "max_pip", "n_mtid", "best_mtid"]
cover.sort(key=lambda r: (r["gene"], r["quant"], -r["n"]))

cov_path = os.path.join(OUTDIR, f"{TRAIT}.eqtl_coverage.tsv")
with open(cov_path, "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in cover:
        fh.write("\t".join(str(r[c]) for c in cols) + "\n")

# best = highest-N dataset that fine-maps a gene under EACH quant method
# (ties broken by max_pip). One row per (gene, quant) with any coverage.
best = {}   # (ensg, quant) -> row
for r in cover:
    k = (r["ensg"], r["quant"])
    cur = best.get(k)
    if cur is None or (r["n"], r["max_pip"]) > (cur["n"], cur["max_pip"]):
        best[k] = r
best_path = os.path.join(OUTDIR, f"{TRAIT}.eqtl_coverage.best.tsv")
with open(best_path, "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in sorted(best.values(), key=lambda r: (r["gene"], r["quant"])):
        fh.write("\t".join(str(r[c]) for c in cols) + "\n")

# ---- 5. console summary ----------------------------------------------------
print("\n== per-gene coverage (quant methods with >=1 fine-mapping dataset) ==")
by_gene = {}
for r in cover:
    by_gene.setdefault(r["ensg"], []).append(r)
covered = 0
for ensg, sym in targets.items():
    rows = by_gene.get(ensg, [])
    if rows:
        covered += 1
        qs = sorted({r["quant"] for r in rows})
        b = sorted(rows, key=lambda r: (-r["n"], -r["max_pip"]))[0]
        print(f"   {sym:7s} {len(rows):3d} dataset-hit(s) via {','.join(qs):20s} "
              f"best: {b['dataset']} {b['study']}/{b['tissue']}({b['quant']}) "
              f"n={b['n']} maxPIP={b['max_pip']}")
    else:
        print(f"   {sym:7s}   0 dataset-hit(s)  (not fine-mapped anywhere)")
print(f"\n[done] {covered}/{len(targets)} genes fine-mapped in \u2265" "1 dataset")
print(f"       -> {cov_path}")
print(f"       -> {best_path}   (one row per gene x quant method -> drives 05a coloc)")
