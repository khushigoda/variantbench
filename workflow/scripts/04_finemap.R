#!/usr/bin/env Rscript
# Step 4 — SuSiE fine-mapping per locus + missense/splice annotation of credible sets.
# snakemake@input$meta / $lead ; snakemake@output$annot
# NOTE: LD is OUT-OF-SAMPLE 1000G EUR (paper used in-sample UKBB LD) -> credible sets can
#       be miscalibrated. DOCUMENT this in the report; it is a scoring point.
suppressMessages({library(data.table); library(susieR)})
set.seed(snakemake@params$seed)   # reproducibility: global seed from config

meta  <- fread(snakemake@input$meta)
lead  <- fread(snakemake@input$lead)
panel <- snakemake@params$panel
mhc   <- snakemake@params$mhc

# ---- TODO ----
# For each locus in lead (skip MHC chr6:mhc$start-mhc$end):
#   1. Subset meta to locus; align to `panel` variants.
#   2. Build LD matrix R from `panel` (plink --r square, or bigsnpr) for locus SNPs.
#   3. z <- beta/se ; fit <- susie_rss(z, R, n=median(n), L=10).
#   4. Extract credible sets + PIP; keep max-PIP variant per CS.
#   5. Annotate CS variants: VEP (Ensembl REST / offline) -> missense/splice consequence;
#      optional SpliceAI. Map coding CS variant -> candidate effector gene.
# Emit annot.tsv: locus, cs, variant_id, pip, gene, consequence.

annot <- data.table(locus=character(), cs=integer(), variant_id=character(),
                    pip=numeric(), gene=character(), consequence=character())
fwrite(annot, snakemake@output$annot, sep = "\t")

# ---- PLOTS (§5 report) ----
# One PIP / regional plot per locus, written INTO the directory() output. The count is
# data-dependent, so we emit a directory rather than named files. Scaffold just creates
# the dir so the rule succeeds; the loop above will fill it with <locus>.pip.png files.
dir.create(snakemake@output$plotdir, showWarnings = FALSE, recursive = TRUE)
