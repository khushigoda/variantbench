#!/usr/bin/env bash
# VariantBench — CAD-only data pull (subset of scripts/download_data.sh).
# TEMPORARY convenience script: pulls ONLY the two CAD outcome files, so you can
# fetch them now without triggering the metabolite / LDSC / eQTL sections.
# Run from the repo root:  bash scripts/download_cad.sh
# Verified live against EBI GWAS Catalog FTP (2026-07).
set -euo pipefail

RAW=data/raw
mkdir -p "$RAW"
base_ftp="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"

# ---- CAD — Aragam et al. 2022 (harmonised GRCh38) ----
#   GCST90132314 = European-only (181,522 cases)  -> ancestry-matched to EUR instruments,
#                  and the accession the authors' released code uses. EUR-only; .315
#                  (multi-ancestry) is NOT downloaded.
# HARMONISED files are already lifted to GRCh38 (genome_assembly: GRCh38 in the
#   -meta.yaml), matching the GRCh38 metabolite exposures -> NO self-liftover needed.
#   Filename pattern: <acc>.h.tsv.gz  (gzipped, ~1.3 GB)
CAD_DIR="$base_ftp/GCST90132001-GCST90133000"
for acc in GCST90132314; do
  echo ">> CAD $acc (harmonised GRCh38)"
  curl -L --fail -C - -o "$RAW/CAD.${acc}.h.tsv.gz" \
     "$CAD_DIR/$acc/harmonised/${acc}.h.tsv.gz"
done

echo "DONE. CAD files in $RAW:"
ls -lh "$RAW"/CAD.GCST90132314.h.tsv.gz
