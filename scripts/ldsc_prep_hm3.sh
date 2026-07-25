#!/usr/bin/env bash
# Pre-restrict a genome-wide meta to HapMap3 SNPs and emit a slim, clean munge input.
# Fixes two failures seen when munging the full genome-wide meta directly:
#   (1) eaf is NA for ~76% of variants (the UKBB-only rows carry no allele frequency).
#       munge does dropna(how="any") across every column it reads, so those variants are
#       discarded BEFORE the HapMap3 merge — a large, needless power loss. LDSC h2/rg do
#       not use allele frequency, so we drop eaf entirely.
#   (2) N (an integer column) is rejected as non-numeric by the LDSC env's old numpy/pandas
#       issubdtype(...,number) check, while the float columns (z) pass. We emit n as a float
#       so it passes the same check. N is genuinely clean (verified across all 103M rows).
# Bonus: one uncompressed ~1.2M-row file -> munge reads it in a single pass, not 42 gzip chunks.
set -euo pipefail
meta="$1"; hm3="$2"; out="$3"
gzip -dc "$meta" | awk -F'\t' -v OFS='\t' '
  FNR==NR { if ($1!="SNP") hm[$1]=1; next }               # load HM3 rsIDs (skip its header)
  FNR==1  { print "SNP","ea","oa","z","n"; next }          # emit slim header
  ($1 in hm) { printf "%s\t%s\t%s\t%s\t%.1f\n", $1,$4,$5,$9,$12 }
' "$hm3" - > "$out"
