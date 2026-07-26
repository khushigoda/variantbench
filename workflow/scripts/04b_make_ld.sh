#!/usr/bin/env bash
# =====================================================================================
# 04b_make_ld.sh  —  Step-4 prep: build per-locus EUR LD matrices for LDL_C fine-mapping
# -------------------------------------------------------------------------------------
# INPUT  : results/loci/LDL_C.lead.tsv          (Step-3 lead variants)
#          results/meta/LDL_C.meta.tsv.gz        (Step-1 meta, for regional GWAS + z)
#          1000G phase3-on-GRCh38 plink2 pgen    (per-chr chr{N}_hg38.*, only needed chrs)
# OUTPUT : data/ref/1kg_phase3_hg38/chr{N}_hg38.{pgen,pvar,psam}  (per-chr panel, cached)
#          results/ld/LDL_C.<gene>.ld              (signed r LD matrix, square)
#          results/ld/LDL_C.<gene>.ld.vars         (variant order for the matrix)
#          results/ld/LDL_C.<gene>.z.tsv           (harmonised, sign-aligned GWAS z)
#          results/ld/LDL_C.ld_manifest.tsv        (per-locus QC: n_var, dropped, flags)
#
# DESIGN (decisions locked with reviewer):
#   * Panel : 1000G phase3-on-GRCh38 HIGH-COV NYGC rebuild (3202-sample) PLINK2 pgen, SPLIT BY CHROMOSOME
#             (chr{N}_hg38 --pfile trios (.pgen + .pvar.zst + .psam)). Provisioned ONCE by
#             scripts/download_data.sh (whole panel, decompressed, "_rs" stripped from pvar
#             names, shared hg38_corrected.psam symlinked per chr). This script only READS
#             the local panel for whichever chromosomes the selected loci touch.
#             Relatedness handled downstream: keep-file restricts to unrelated EUR FOUNDERS
#             (525 of 633 EUR; 108 trio offspring dropped via PAT/MAT flags) so LD is unbiased.
#   * Pop   : EUR only. Keep file built by INSPECTING the .psam header (ID + SuperPop
#             columns detected at runtime — schema differs across 1000G releases).
#   * Loci  : eligibility screen on the 321 leads, THEN curated CAD-effector genes.
#             Screen = drop maf_class=="unknown" (0 now), drop MHC chr6:25-34Mb (6),
#             drop rare leads MAF<0.1% (14 — a MAF>0.1% panel can't support them).
#             Curation = nearest strongest eligible lead within 500 kb of each effector
#             gene anchor (gene->coord table below). One lead per gene.
#   * Window: FIXED +/-500 kb for every locus (no auto-widen). Long-range-LD truncation
#             (e.g. APOE) is a documented Step-4 MCQ point, flagged in the manifest.
#   * Harmon: match panel<->GWAS on chr:pos + allele pair; allow ref/alt SWAP (flip z
#             sign); DROP strand-ambiguous palindromes (A/T, C/G). LD matrix rows and
#             the z-vector end in identical position order => same variants, same allele.
#   * LD     : SIGNED r (not r^2) — SuSiE needs the signed correlation matrix.
#
# RUN    : this script is normally driven by the Snakemake rule `make_ld` (trait = wildcard,
#          all knobs = config). It also runs standalone via env vars (defaults in [brackets]):
#            conda activate variantbench   (provides bcftools, awk, Rscript)
#            plink2 comes from its own Rosetta/x86 env: PLINK2="conda run -n plink2 plink2"
#            TRAIT=LDL_C bash workflow/scripts/04b_make_ld.sh          # curated CAD loci [default]
#            TRAIT=LDL_C CURATE=0 bash .../04b_make_ld.sh              # every eligible lead
#            TRAIT=LDL_C LOCI=PCSK9 bash .../04b_make_ld.sh            # single-locus sample test
#            TRAIT=LDL_C ANCHORS=config/other_anchors.tsv bash ...     # different curation set
# =====================================================================================
set -euo pipefail

# ---------------------------------------------------------------- 0. config / paths ---
# All knobs come from the environment (Snakemake sets them from config+wildcards); each has
# a standalone default. TRAIT is the only required one — nothing else is trait-hardcoded.
TRAIT="${TRAIT:?set TRAIT=<metabolite> (e.g. LDL_C) — nothing is hardcoded}"
LEAD_TSV="results/loci/${TRAIT}.lead.tsv"
META_GZ="results/meta/${TRAIT}.meta.tsv.gz"
REFDIR="${REFDIR:-data/ref/1kg_phase3_hg38}"     # per-chromosome panel cache (shared across traits)
OUTDIR="${OUTDIR:-results/ld}"
MANIFEST="${OUTDIR}/${TRAIT}.ld_manifest.tsv"
LOCI_TSV="${OUTDIR}/${TRAIT}.loci.tsv"

WINDOW="${WINDOW:-500000}"              # +/- fine-mapping window (bp)
SELECT_RADIUS="${SELECT_RADIUS:-500000}"  # a gene anchor claims the strongest eligible lead within this
MAF_MIN="${MAF_MIN:-0.001}"            # panel QC floor (MAF > 0.1%) — matches Step-3 common-variant floor
EDGE_FLAG_BP="${EDGE_FLAG_BP:-50000}"  # manifest flag: signal within this many bp of window edge

# Curation switch: CURATE=1 [default] -> keep the strongest eligible lead near each gene anchor
# (ANCHORS file). CURATE=0 -> every eligible lead becomes its own locus (future: whole lead list).
CURATE="${CURATE:-1}"
ANCHORS="${ANCHORS:-config/cad_effector_anchors.tsv}"   # gene<TAB>chr<TAB>anchor_bp ; swap per trait

# R half lives in a separate script (04a_ld_helpers.R) so it can run under the R env while
# this orchestrator + plink2 run under the plink env. RSCRIPT is the Rscript to use (system
# Rscript by default — confirmed on PATH in every env; override to the R env's if preferred).
RSCRIPT="${RSCRIPT:-Rscript}"
RHELP="${RHELP:-workflow/scripts/04a_ld_helpers.R}"
# plink2 lives in its own Rosetta/x86 conda env on arm64 Macs (no osx-arm64 bioconda build).
# Default: reach it via `conda run -n plink2 plink2`; override PLINK2 if it's on PATH elsewhere.
PLINK2="${PLINK2:-conda run -n plink2 plink2}"
[[ -f "$RHELP" ]] || { echo "[STOP] R helper '$RHELP' not found." >&2; exit 5; }

# The 1000G phase3-GRCh38 panel (chr{N}_hg38 --pfile trios (.pgen + .pvar.zst + .psam)) is provisioned
# ONCE by scripts/download_data.sh into $REFDIR. This script only READS it — no downloads here.

mkdir -p "$OUTDIR"

echo "=========================================================================="
echo " 04b_make_ld.sh — Step-4 LD prep for ${TRAIT}"
echo " window=+/-${WINDOW}bp  maf_min=${MAF_MIN}  panel=1000G-phase3-GRCh38(525 unrel EUR founders)"
echo "=========================================================================="

# ================================================ 1. locus table (screen + curate) =====
# BUILT FIRST: the selected loci determine which chromosomes we need from the local panel.
# Eligibility screen + CAD-effector curation, entirely from the lead table.
echo "[1/5] selecting loci (eligibility screen; curate=${CURATE}) ..."
[[ "$CURATE" == "1" && ! -f "$ANCHORS" ]] && { echo "[STOP] CURATE=1 but anchor file '$ANCHORS' not found." >&2; exit 4; }
# R half runs under the R env's Rscript (base R only). LOCI (optional comma-list) restricts curation.
"$RSCRIPT" "$RHELP" select "$LEAD_TSV" "$LOCI_TSV" "$WINDOW" "$SELECT_RADIUS" "$CURATE" "$ANCHORS" "${LOCI:-}"

# ============================ 2. verify the local panel has the chromosomes we need =====
# The locus table names the chromosomes; the panel was provisioned by download_data.sh.
CHRS=$(tail -n +2 "$LOCI_TSV" | cut -f2 | sort -un | tr '\n' ' ')
echo "[2/5] chromosomes required by selected loci: ${CHRS}"
[[ -f "${REFDIR}/hg38_corrected.psam" ]] || { echo "[STOP] ${REFDIR}/hg38_corrected.psam missing — run scripts/download_data.sh first." >&2; exit 3; }
for chr in $CHRS; do
  base="${REFDIR}/chr${chr}_hg38"
  # trio is: .pgen (decompressed) + .pvar.zst (kept compressed, read via vzs) + .psam
  for ext in pgen pvar.zst psam; do
    [[ -e "${base}.${ext}" ]] || { echo "[STOP] panel file ${base}.${ext} missing — run scripts/download_data.sh first." >&2; exit 3; }
  done
done
echo "[2/5] local panel present for all required chromosomes."

# ================================================= 3. EUR keep file (inspect .psam) ====
# .psam schema varies (#FID+IID vs #IID; SuperPop column name). Detect at runtime.
KEEP="${REFDIR}/1kg_eur.keep"
echo "[3/5] building EUR keep file from hg38_corrected.psam header ..."
"$RSCRIPT" "$RHELP" keep "${REFDIR}/hg38_corrected.psam" "$KEEP"

# ============================ 4-5. per-locus: extract -> harmonise -> signed LD matrix =
# EUR filter (--keep) is applied per-locus at extraction — no whole-panel pre-subset needed.
echo "[4-5/5] per-locus extract (EUR keep) + harmonise + signed LD ..."
: > "$MANIFEST"
echo -e "gene\tchr\tlead_snp\tlead_pos\twindow_bp\tn_panel\tn_gwas\tn_matched\tn_swapped\tn_palindrome_dropped\tn_final\tedge_flag\tstatus" >> "$MANIFEST"

# iterate loci (skip header)
tail -n +2 "$LOCI_TSV" | while IFS=$'\t' read -r gene chr lead_snp lead_pos ea oa z nlp start end; do
  tag="${TRAIT}.${gene}"
  echo "    -- ${gene}  chr${chr}:${start}-${end}  (lead ${lead_snp}) --"

  # 4. extract locus from per-chr panel, EUR-only: biallelic SNPs, MAF>0.1%, region-clipped
  #    vzs: pvar is stored compressed (chr{N}_hg38.pvar.zst); plink2 reads it in place.
  #    --rm-dup exclude-all: the cog-genomics build split multiallelics into biallelic rows that
  #    keep the SAME rsID, so --snps-only/--max-alleles 2 still lets duplicate-ID rows through.
  #    A split-multiallelic rsID is ambiguous (can't map the GWAS z to a single allele pair), so
  #    drop ALL copies — this keeps the LD matrix and the z-vector on one unambiguous variant set.
  $PLINK2 --pfile "${REFDIR}/chr${chr}_hg38" vzs --keep "$KEEP" \
         --chr "$chr" --from-bp "$start" --to-bp "$end" \
         --snps-only --max-alleles 2 --maf "$MAF_MIN" --rm-dup exclude-all \
         --make-pgen --out "${OUTDIR}/${tag}.panel" >/dev/null 2>&1 || {
           echo -e "${gene}\t${chr}\t${lead_snp}\t${lead_pos}\t$((2*WINDOW))\t0\t0\t0\t0\t0\t0\tNA\tEMPTY_PANEL" >> "$MANIFEST"; continue; }

  # 5. regional GWAS slice (meta: chr=$2 pos=$3 ea=$4 oa=$5 beta=$7 se=$8 z=$9 nlp=$11)
  zcat < "$META_GZ" | awk -F'\t' -v c="$chr" -v s="$start" -v e="$end" \
    'NR==1{print "chr\tpos\tea\toa\tz\tnlp"; next} $2==c && $3>=s && $3<=e {print $2"\t"$3"\t"$4"\t"$5"\t"$9"\t"$11}' \
    > "${OUTDIR}/${tag}.gwas.tsv"

  # 5b. harmonise panel<->GWAS -> aligned z + extract list (position-ordered) — R half
  "$RSCRIPT" "$RHELP" harmonize "${OUTDIR}/${tag}.panel.pvar" "${OUTDIR}/${tag}.gwas.tsv" \
             "${OUTDIR}/${tag}" "$gene" "$chr" "$lead_pos" "$WINDOW" "$EDGE_FLAG_BP" "$MANIFEST"

  # fill lead_snp into the manifest row just written (R wrote it blank in col3)
  # (cosmetic; the gene label is the join key downstream)

  # 6. signed LD matrix on the harmonised set (r, NOT r^2), square
  if [[ -s "${OUTDIR}/${tag}.extract.ids" ]] && [[ $(wc -l < "${OUTDIR}/${tag}.extract.ids") -ge 2 ]]; then
    $PLINK2 --pfile "${OUTDIR}/${tag}.panel" --extract "${OUTDIR}/${tag}.extract.ids" \
           --r-unphased square ref-based \
           --out "${OUTDIR}/${tag}" >/dev/null 2>&1 || echo "       [warn] plink2 --r failed for ${gene}"
    # plink2 writes <tag>.unphased.vcor1 (matrix) + <tag>.unphased.vcor1.vars (order).
    # Normalise names so Step-4 finds them predictably:
    [[ -f "${OUTDIR}/${tag}.unphased.vcor1" ]]      && mv "${OUTDIR}/${tag}.unphased.vcor1"      "${OUTDIR}/${tag}.ld"
    [[ -f "${OUTDIR}/${tag}.unphased.vcor1.vars" ]] && mv "${OUTDIR}/${tag}.unphased.vcor1.vars" "${OUTDIR}/${tag}.ld.vars"
    # align z-vector to the LD matrix's own variant order + assert the ID sets match exactly
    # (SuSiE needs z[i] <-> R[i,i]; guarantee it by construction rather than by coincident sorting)
    "$RSCRIPT" "$RHELP" alignz "${OUTDIR}/${tag}.ld.vars" "${OUTDIR}/${tag}.z.tsv" || {
      echo "       [FAIL] z/LD alignment mismatch for ${gene} — see error above"; exit 1; }
    echo "       wrote ${tag}.ld  (+ .ld.vars, .z.tsv, z aligned to LD order)"
  else
    echo "       [skip LD] too few harmonised variants for ${gene}"
  fi

  # tidy per-locus scratch (keep .ld, .ld.vars, .z.tsv, .gwas.tsv)
  rm -f "${OUTDIR}/${tag}.panel."{pgen,pvar,psam,log} 2>/dev/null || true
done

# ================================================================== done / manifest =
echo "[done] LD prep complete. Manifest:"
column -t -s$'\t' "$MANIFEST" 2>/dev/null || cat "$MANIFEST"
echo
echo "Per-locus SuSiE inputs (Step 4): results/ld/${TRAIT}.<gene>.{ld,ld.vars,z.tsv}"
echo "LD matrix = signed r (525-EUR-founder out-of-sample reference); z = meta, allele-aligned to panel ALT."
