#!/usr/bin/env bash
# Step 0 (normalize/prepare) — split a raw GWAS-SSF .tsv.gz into per-chromosome
# files in ONE decompression pass, selecting only the columns Step 1 needs.
#
# Why: an 8 GB machine cannot hold two 26M/96M-row cohort tables + their merge.
# gzip is not seekable, so a per-chromosome loop that re-reads the whole .tsv.gz
# per chromosome would decompress it 23x. Instead we decompress ONCE and route
# each row to its chromosome's gzip stream (awk keeps 23 pipes open in parallel).
# Step 1 then reads one small chr-file per cohort at a time -> peak memory = 1 chr.
#
# Usage: 00_split_by_chr.sh <raw.tsv.gz> <outdir> <comma-sep 1-based col idxs>
#   col 1 (chromosome) MUST be included and is used for routing.
set -euo pipefail
IN="$1"; OUT="$2"; COLS="$3"
mkdir -p "$OUT"
# clear any partial previous run
rm -f "$OUT"/chr*.tsv.gz "$OUT"/.done 2>/dev/null || true

gzip -dc "$IN" | awk -v cols="$COLS" -v outdir="$OUT" '
BEGIN{ n=split(cols, C, ",") }
{
  line=$(C[1]); for(k=2;k<=n;k++) line=line "\t" $(C[k])
  if(NR==1){ header=line; next }
  c=$1
  f=outdir"/chr"c".tsv.gz"
  cmd="gzip -c > \"" f "\""
  if(!(c in seen)){ print header | cmd; seen[c]=1 }
  print line | cmd
}'
touch "$OUT/.done"
echo "split $(basename "$IN") -> $OUT/  ($(ls "$OUT"/chr*.tsv.gz 2>/dev/null | wc -l | tr -d ' ') chr files)"
