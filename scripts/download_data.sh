#!/usr/bin/env bash
# VariantBench — data acquisition
# All accessions/paths below were verified live against the GWAS Catalog REST API
# and EBI FTP (2026-07). Metabolites chosen: Total BCAA + Lactate + glucose + LDL_C.
# Cohort decoding (from GWAS Catalog sample sizes):
#   EstBB    = 185,352 EUR
#   UKBB_EUR = 413,897 EUR
#   meta_EUR = 599,249 EUR   (= EstBB + UKBB_EUR; the paper's target meta-analysis)
set -euo pipefail

RAW=data/raw
REF=data/ref
mkdir -p "$RAW" "$REF"
base_ftp="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"

# ---- section selector -------------------------------------------------------
# Run only the sections you need (default = all). Sections: metab cad ldsc panel eqtl.
#   bash scripts/download_data.sh panel        # only the LD panel (Step 4)
#   bash scripts/download_data.sh metab ldsc   # exposures + LDSC reference
#   bash scripts/download_data.sh              # everything
# Guards mean re-running never re-fetches a complete file, but a FULL run still re-verifies
# the multi-GB sumstats (gzip -t, several minutes); selecting `panel` skips all of that.
SECTIONS="${*:-all}"
want () { [[ " $SECTIONS " == *" all "* || " $SECTIONS " == *" $1 "* ]]; }
# -----------------------------------------------------------------------------

# ---- completeness checks ----------------------------------------------------
# The old guard was `[[ -s FILE ]]` (non-empty). That treats a truncated/interrupted
# download as "present" and skips it forever. These helpers verify a file is a COMPLETE
# download before we skip it; anything that fails is re-fetched (curl -C - resumes).
#
# gzip_ok: the whole gzip stream decompresses cleanly (catches truncation at the byte
#   level — a half-written .gz fails here). Cheap: decompresses to /dev/null, no temp file.
gzip_ok () { [[ -s "$1" ]] && gzip -t "$1" 2>/dev/null; }

# md5_ok: if a sidecar <base>.meta.yaml carries a data_file_md5sum, verify against it.
#   Returns success (0) when there is NO checksum to check, so it never blocks on absence.
md5_ok () {  # $1 = data file, $2 = meta.yaml (optional)
  local f="$1" meta="${2:-}"
  [[ -f "$meta" ]] || return 0
  local want; want=$(grep -Eo 'data_file_md5sum:[[:space:]]*[0-9a-fA-F]{32}' "$meta" \
                     | grep -Eo '[0-9a-fA-F]{32}' | head -1)
  [[ -n "$want" ]] || return 0
  local got; got=$( (md5sum "$f" 2>/dev/null || md5 -q "$f") | awk '{print $1}')
  [[ "$got" == "$want" ]]
}

# complete_gz: skip-if-true guard for a gzipped data file (+ optional meta for md5).
complete_gz () { gzip_ok "$1" && md5_ok "$1" "${2:-}"; }
# -----------------------------------------------------------------------------

ftp_url () {  # $1 = accession -> full FTP file URL
  local acc="$1"
  local n=${acc:4}
  local lo=$(( ((10#$n - 1)/1000)*1000 + 1 ))
  local hi=$(( lo + 999 ))
  printf "%s/GCST%08d-GCST%08d/%s/%s.tsv.gz" \
         "$base_ftp" "$lo" "$hi" "$acc" "$acc"
}

# ---- 1. Metabolite exposure sumstats (Tambets et al., PMID 42162431) ----
# EUR-only, 4 traits x 3 cohorts. Accessions from Supplementary Table S12
# ("Download paths for GWAS summary statistics hosted at the GWAS Catalog").
# Cohort scheme (EUR only): EstBB, UKBB_EUR, metaEUR (= EstBB+UKBB_EUR).
# The non-EUR UKBB ancestries and meta_full are NOT pulled.
#   reported_trait in S12:
#     BCAA    = "Total concentration of branched-chain amino acids (leucine + isoleucine + valine)"
#     LDL_C   = "LDL cholesterol"   (NOT "Clinical LDL cholesterol", a separate trait)
#     Lactate = "Lactate"
#     Glucose = "Glucose"
# name:accession
METAB=(
  "BCAA_EstBB:GCST90449544"
  "BCAA_UKBB_EUR:GCST90449793"
  "BCAA_metaEUR:GCST90451287"
  "Lactate_EstBB:GCST90449404"
  "Lactate_UKBB_EUR:GCST90449653"
  "Lactate_metaEUR:GCST90451147"
  "LDL_C_EstBB:GCST90449408"
  "LDL_C_UKBB_EUR:GCST90449657"
  "LDL_C_metaEUR:GCST90451151"
  "Glucose_EstBB:GCST90449379"
  "Glucose_UKBB_EUR:GCST90449628"
  "Glucose_metaEUR:GCST90451122"
)
if want metab; then
for entry in "${METAB[@]}"; do
  name="${entry%%:*}"; acc="${entry##*:}"
  gz="$RAW/${name}.${acc}.tsv.gz"; meta="$RAW/${name}.${acc}.meta.yaml"
  # fetch metadata first so md5 verification is possible on this run
  if [[ ! -f "$meta" ]]; then
    curl -L --fail -o "$meta" "$(ftp_url "$acc")-meta.yaml" \
      || echo "WARNING: metadata unavailable for $acc" >&2
  fi
  if complete_gz "$gz" "$meta"; then
    echo ">> $name ($acc): already present & complete, skipping"
    continue
  fi
  [[ -s "$gz" ]] && echo ">> $name ($acc): incomplete/failed check — re-fetching (resume)"
  url=$(ftp_url "$acc")
  echo ">> $name  ($acc)"
  curl -L --fail -C - -o "$gz" "$url"
  complete_gz "$gz" "$meta" || { echo "ERROR: $name ($acc) still fails completeness check" >&2; exit 1; }
done
else echo ">> [skip] metabolite exposures (section not selected)"; fi

# ---- 2. Disease OUTCOME sumstats for MR ----
# CAD — Aragam et al. 2022. EUR-only:
#   GCST90132314 = European-only (181,522 cases)  -> ancestry-matched to EUR instruments,
#                  and what the authors' released code actually used.
#   (GCST90132315, the multi-ancestry accession the paper's Methods cite, is NOT downloaded
#    — EUR-only is ancestry-matched to our EUR instruments and is the primary/only outcome.)
# We download the GWAS Catalog HARMONISED file, which is already lifted to GRCh38
#   (genome_assembly: GRCh38 in the -meta.yaml). Metabolite exposures are GRCh38, so
#   pulling harmonised CAD means NO self-liftover step is needed. This is also what the
#   authors effectively did from their code (their Aragam_2022_..._harmonized.parquet).
#   Filename pattern: <acc>.h.tsv.gz  (gzipped, ~1.3 GB)
if want cad; then
CAD_DIR="$base_ftp/GCST90132001-GCST90133000"
for acc in GCST90132314; do
  gz="$RAW/CAD.${acc}.h.tsv.gz"
  if gzip_ok "$gz"; then
    echo ">> CAD $acc: already present & complete, skipping"
    continue
  fi
  [[ -s "$gz" ]] && echo ">> CAD $acc: incomplete — re-fetching (resume)"
  echo ">> CAD $acc (harmonised GRCh38)"
  curl -L --fail -C - -o "$gz" \
     "$CAD_DIR/$acc/harmonised/${acc}.h.tsv.gz"
  gzip_ok "$gz" || { echo "ERROR: CAD $acc still fails gzip check" >&2; exit 1; }
done

# T2D — Suzuki et al. 2024 (DIAMANTE / T2DGGI). MANUAL download (agreement required),
# NOT scriptable via curl. EUR-only, from the DIAGRAM downloads page:
#   https://www.diagram-consortium.org/downloads.html
#   - EUR_Metal_LDSC-CORR_Neff.v2.txt  (European-only)   -> primary (ancestry-matched)
#   (All_Metal_LDSC-CORR_Neff.v2.txt, multi-ancestry, is NOT downloaded — EUR-only only.)
# BUILD NOTE: DIAGRAM serves these ONLY in GRCh37/hg19 (confirmed in their
#   T2DGGI_GWAS_summary_statistics.pdf: "chromosome and position (hg19, build 37)").
#   There is NO pre-harmonised GRCh38 download for T2D (unlike CAD on the GWAS Catalog).
#   => T2D STILL requires GRCh37 -> GRCh38 liftover in Step 0 before downstream analysis.
#   The authors did their own liftover here it seems from code (their Suzuki_2024_..._harmonized.parquet).
# HEADER NOTE: first column is spelled "Chromsome" (sic) in the source files;
#   00_normalize.R must map that exact string.
#
# After downloading, place + rename + compress (run from repo root):
#   mv ~/Downloads/EUR_Metal_LDSC-CORR_Neff.v2.txt $RAW/T2D.Suzuki2024.EUR.tsv
#   gzip $RAW/T2D.Suzuki2024.EUR.tsv   # -> T2D.Suzuki2024.EUR.tsv.gz
echo ">> T2D: place manually-downloaded DIAGRAM file (GRCh37):"
echo "        $RAW/T2D.Suzuki2024.EUR.tsv.gz (EUR-only, primary)"
# manual file — verify its gzip integrity if present, warn (don't fail) if absent
if [[ -e "$RAW/T2D.Suzuki2024.EUR.tsv.gz" ]]; then
  gzip_ok "$RAW/T2D.Suzuki2024.EUR.tsv.gz" \
    && echo "   T2D EUR file: present & gzip-complete" \
    || echo "   WARNING: T2D EUR file present but FAILS gzip check — re-download manually" >&2
else
  echo "   (T2D EUR file not yet placed — manual step)"
fi
else echo ">> [skip] CAD/T2D outcome sumstats (section not selected)"; fi

# ---- 3. LDSC reference data (Step 2) ----
# HapMap3 SNP list (for munge --merge-alleles) and EUR LD scores (for --ref-ld/--w-ld).
# Both guards verify a COMPLETE download (gzip stream valid), then materialize the
# uncompressed forms the config points at: data/ref/w_hm3.snplist (file) and
# data/ref/eur_w_ld_chr/ (DIRECTORY — the tar must be extracted, not just downloaded).
if want ldsc; then
if ! gzip_ok "$REF/w_hm3.snplist.gz"; then
  [[ -s "$REF/w_hm3.snplist.gz" ]] && echo ">> w_hm3.snplist.gz: incomplete — re-fetching (resume)"
  curl -L --fail -C - -o "$REF/w_hm3.snplist.gz" \
    "https://zenodo.org/records/7773502/files/w_hm3.snplist.gz?download=1"
  gzip_ok "$REF/w_hm3.snplist.gz" || { echo "ERROR: w_hm3.snplist.gz fails gzip check" >&2; exit 1; }
fi

if ! gzip_ok "$REF/eur_w_ld_chr.tar.gz"; then
  [[ -s "$REF/eur_w_ld_chr.tar.gz" ]] && echo ">> eur_w_ld_chr.tar.gz: incomplete — re-fetching (resume)"
  curl -L --fail -C - -o "$REF/eur_w_ld_chr.tar.gz" \
    "https://zenodo.org/records/8182036/files/eur_w_ld_chr.tar.gz?download=1"
  gzip_ok "$REF/eur_w_ld_chr.tar.gz" || { echo "ERROR: eur_w_ld_chr.tar.gz fails gzip check" >&2; exit 1; }
fi

# uncompressed HapMap3 list (config: hm3_snplist)
if [[ ! -s "$REF/w_hm3.snplist" ]]; then
  gunzip -k "$REF/w_hm3.snplist.gz"
fi

# extracted LD-score directory (config: ld_scores = data/ref/eur_w_ld_chr/).
# LDSC needs the 22 per-chr .l2.ldscore.gz + .M_5_50 files, so extract if the dir is
# missing or looks incomplete (<22 ldscore files).
if [[ ! -d "$REF/eur_w_ld_chr" ]] || \
   [[ $(find "$REF/eur_w_ld_chr" -name '*.l2.ldscore.gz' 2>/dev/null | wc -l) -lt 22 ]]; then
  echo ">> extracting eur_w_ld_chr.tar.gz"
  tar -xzf "$REF/eur_w_ld_chr.tar.gz" -C "$REF"
fi
else echo ">> [skip] LDSC reference data (section not selected)"; fi

# ---- 4. LD reference panel — 1000G phase3 GRCh38 pgen (Step 4 fine-mapping) ----
# 1000 Genomes phase 3, GRCh38, split-by-chromosome, rsID-annotated (dbSNP156),
# KING-pedigree-corrected. Source: cog-genomics.org/plink/2.0/resources.
# The 45 versioned Dropbox URLs (22 pgen + 22 _rs.pvar + 1 shared corrected .psam) live
# in config/panel_urls.tsv (extracted from the resources page; they change over time —
# re-extract if a link 404s). We pull the WHOLE panel here so the Step-4 analysis script
# is pure "read local data + analyze" with no download logic. ~5.6 GB compressed.
#
# PLINK2 .zst handling (per the resource page's own note):
#   * .pgen.zst  MUST be decompressed before use  -> `plink2 --zst-decompress` (also our
#                integrity gate: a truncated .zst fails here, so we never skip a bad file).
#   * .pvar.zst  does NOT need decompressing  -> `--pfile <base> vzs` reads it compressed.
#                We keep it zipped (saves ~4 GB across 22 chr) and just give it the base
#                name PLINK expects: chr{N}_hg38.pvar.zst (dropping the "_rs" tag).
# We deliberately do NOT depend on a standalone `zstd` binary (may be absent) — plink2 is
# the only decompressor used. Each chr ends up as a `--pfile chr{N}_hg38 vzs` trio:
# chr{N}_hg38.pgen (plain) + chr{N}_hg38.pvar.zst (compressed) + chr{N}_hg38.psam (symlink).
# CHROMOSOME SCOPE: we only need the chromosomes our fine-mapping loci touch, not all 22.
#   PANEL_CHRS defaults to the unique chromosomes in the anchor file (config/cad_effector_anchors.tsv)
#   — the curated CAD set spans chr 1,2,5,6,7,8,9,11,16,19 (10 of 22, ~2.9 GB not ~5.6 GB).
#   Override for a wider run:  PANEL_CHRS="1 2 3 ... 22"  or  PANEL_CHRS=all
# Only the pgen/rs-pvar for the selected chromosomes (plus the one shared psam) are fetched;
# every other line in panel_urls.tsv is skipped.
PANEL_DIR="$REF/1kg_phase3_hg38"
PANEL_TSV="config/panel_urls.tsv"
ANCHORS_TSV="${ANCHORS:-config/cad_effector_anchors.tsv}"
if want panel; then
mkdir -p "$PANEL_DIR"
if [[ ! -f "$PANEL_TSV" ]]; then
  echo "WARNING: $PANEL_TSV missing — skipping LD panel download (Step 4 needs it)." >&2
else
  command -v plink2 >/dev/null 2>&1 || { echo "ERROR: plink2 not on PATH (conda activate variantbench)"; exit 1; }

  # resolve which chromosomes to pull
  if [[ "${PANEL_CHRS:-}" == "all" ]]; then
    CHRS=$(seq 1 22)
  elif [[ -n "${PANEL_CHRS:-}" ]]; then
    CHRS="$PANEL_CHRS"
  elif [[ -f "$ANCHORS_TSV" ]]; then
    CHRS=$(grep -v '^#' "$ANCHORS_TSV" | awk -F'\t' 'NF>=3{print $2}' | sort -n -u | tr '\n' ' ')
  else
    echo "WARNING: no PANEL_CHRS and no anchor file '$ANCHORS_TSV' — defaulting to all 22." >&2
    CHRS=$(seq 1 22)
  fi
  echo ">> LD panel: fetching chromosomes -> $CHRS"

  # build the set of filenames we actually want (pgen + rs-pvar for CHRS, plus the shared psam)
  # portable (bash 3.2 / macOS default): a space-delimited string, no associative arrays.
  WANT_FILE=" "
  for c in $CHRS; do
    WANT_FILE="${WANT_FILE}chr${c}_hg38.pgen.zst chr${c}_hg38_rs.pvar.zst "
  done
  # fetch only the wanted files (curl --fail rejects HTML error pages; -C - resumes partials)
  while IFS=$'\t' read -r fname url; do
    [[ "$fname" =~ ^#|^$ ]] && continue
    # always take the shared psam; otherwise only files in the wanted set
    if [[ "$fname" != *.psam && "$WANT_FILE" != *" $fname "* ]]; then continue; fi
    out="$PANEL_DIR/$fname"
    [[ -s "$out" ]] && { echo ">> panel $fname: present, skipping download"; continue; }
    echo ">> panel $fname"; curl -L --fail -C - -o "$out" "$url" \
      || { echo "ERROR: download failed for $fname (link may be stale in $PANEL_TSV)"; exit 1; }
  done < "$PANEL_TSV"

  # normalize each selected chromosome into a `--pfile <base> vzs` trio
  for c in $CHRS; do
    base="$PANEL_DIR/chr${c}_hg38"
    # pgen: MUST decompress (also validates the download; re-run resumes a partial then passes)
    [[ -f "${base}.pgen" ]] || plink2 --zst-decompress "${base}.pgen.zst" "${base}.pgen" \
      || { echo "ERROR: ${base}.pgen.zst failed to decompress (truncated? re-run to resume)"; exit 1; }
    # pvar: KEEP COMPRESSED, just give it the base name (--pfile ... vzs reads .pvar.zst)
    [[ -e "${base}.pvar.zst" ]] || ln -sf "chr${c}_hg38_rs.pvar.zst" "${base}.pvar.zst"
    # psam: one shared corrected file for every chromosome
    [[ -e "${base}.psam" ]]     || ln -sf "hg38_corrected.psam"       "${base}.psam"
  done
  echo ">> LD panel ready for chr [$CHRS]: $PANEL_DIR/chr{N}_hg38.{pgen,pvar.zst,psam}  (read with --pfile ... vzs)"
fi
else echo ">> [skip] LD reference panel (section not selected)"; fi

# ---- 5. Molecular QTLs for colocalization (Step 5) ----
# eQTL Catalogue (fine-mapped SuSiE credible sets / LBFs) for the coloc locus.
# Pick the tissue/dataset once you know the locus (e.g. blood eQTL for GP6/lactate
# or liver/skeletal-muscle for BCAA pathway). Browse:
#   https://www.ebi.ac.uk/eqtl/  (FTP: https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/)
if want eqtl; then
echo ">> eQTL Catalogue: select dataset for your coloc locus (see data/README.md)"
fi

echo "DONE (sections: $SECTIONS). Verify md5sums against each *.meta.yaml (data_file_md5sum field)."
