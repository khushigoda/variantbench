#!/usr/bin/env bash
# Bridge a GWAS to a munge-ready TSV (gzipped), PRE-RESTRICTED to HapMap3.
#
# WHY pre-restrict here instead of letting munge do it from the full file: munge_sumstats.py
# (LDSC's frozen Python-2.7 pandas) parses EVERY row into a DataFrame before --merge-alleles
# discards the ~99% that aren't HapMap3. On the genome-wide cohorts (26-104M variants) that
# parse alone runs ~40 min to several HOURS per file on an 8 GB box, silently (munge prints
# nothing until it finishes). Filtering to the ~1.2M HapMap3 rsIDs FIRST — a cheap streaming
# awk join against w_hm3.snplist — means munge reads a ~1.2M-row file in one fast pass
# (~1-2 min) and prints its banner immediately. The scientific result is IDENTICAL: this
# keeps only rsIDs that would survive --merge-alleles anyway, and --merge-alleles STILL runs
# (allele harmonisation + strand-ambiguous/mismatch drops), as do the FRQ/INFO QC defaults.
# We just move the HapMap3 row-count cut one step earlier, off munge's slow parser.
#
# This is otherwise a pure FORMAT bridge (no MAF/INFO filtering — munge's own defaults do that).
# It also fixes the three format incompatibilities munge cannot handle itself:
#   1. P is stored as neg_log10(P) in every source file -> emit a real P. A 1e-300 floor
#      keeps ultra-significant loci from underflowing to 0 (which would make Z = Inf when
#      munge derives Z = sign(BETA) * |qnorm(P/2)|).
#   2. N must be a FLOAT. The LDSC env's Python-2.7 pandas infers an integer column as
#      object dtype and rejects it ("Columns ['N'] are expected to be numeric").
#   3. Column names -> canonical SNP/A1/A2/BETA/P/N/FRQ/INFO so munge auto-detects them.
#      Naming the freq column FRQ triggers the MAF default; naming it INFO triggers the
#      INFO default. (We still pass explicit --snp/--a1/--a2/--N-col flags to be safe.)
#
# Signed statistic = BETA (matches the authors: they pass BETA and let munge pair it with P).
#
#   kind=raw  : per-cohort GWAS-SSF (EstBB / UKBB_EUR). 1-based cols:
#               effect_allele=3 other_allele=4 beta=5 eaf=7 neg_log10p=8 rsid=9 info=11 n=13
#               -> SNP A1 A2 BETA P N FRQ INFO   (INFO present -> --info-min 0.9 applies)
#   kind=meta : our Step-1 meta (results/meta/<m>.meta.tsv.gz). 1-based cols:
#               SNP(rsid)=1 ea=4 oa=5 eaf=6 beta=7 neg_log10p=11 n=12
#               -> SNP A1 A2 BETA P N FRQ        (no INFO col -> INFO filter simply n/a;
#                                                 HM3 + MAF only, as agreed for the meta)
set -euo pipefail
kind="$1"; src="$2"; out="$3"; hm3="$4"   # hm3 = w_hm3.snplist (3-col: SNP A1 A2, rsID-keyed)

# awk reads TWO inputs: (1) the HapMap3 snplist -> build an rsID set; (2) the decompressed
# GWAS stream -> keep only rows whose SNP (rsID) is in that set, reformat, convert P, float N.
awk -F'\t' -v OFS='\t' -v kind="$kind" '
  # neg_log10(P) -> P, floored at 1e-300 so top loci never underflow to exactly 0.
  function p_from_nlog(x,   p) {
    if (x == "NA" || x == "") return "NA"
    if (x + 0 > 300) return "1e-300"
    p = exp(-(x + 0) * log(10))
    if (p < 1e-300) return "1e-300"
    return p
  }
  # First input = HapMap3 snplist: load col1 (rsID) into hm3[], skip its header.
  FNR == NR { if (FNR > 1) hm3[$1] = 1; next }
  # Second input = GWAS. Row 1 is the GWAS header: emit our output header instead.
  FNR == 1 {
    if (kind == "raw") print "SNP","A1","A2","BETA","P","N","FRQ","INFO"
    else               print "SNP","A1","A2","BETA","P","N","FRQ"
    next
  }
  {
    if (kind == "raw") { snp=$9; a1=$3; a2=$4; b=$5; nl=$8; n=$13; frq=$7; info=$11 }
    else               { snp=$1; a1=$4; a2=$5; b=$7; nl=$11; n=$12; frq=$6 }
    if (!(snp in hm3)) next                                   # HapMap3 restriction
    if (snp=="NA" || snp=="" || b=="NA" || b=="" || n=="NA" || n=="") next
    p = p_from_nlog(nl); if (p=="NA") next
    if (kind == "raw")
      printf "%s\t%s\t%s\t%s\t%s\t%.1f\t%s\t%s\n", snp, toupper(a1), toupper(a2), b, p, n, frq, info
    else
      printf "%s\t%s\t%s\t%s\t%s\t%.1f\t%s\n",     snp, toupper(a1), toupper(a2), b, p, n, frq
  }
' "$hm3" <(gzip -dc "$src") | gzip > "$out"
