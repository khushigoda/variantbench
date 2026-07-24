#!/usr/bin/env Rscript
# Step 7 — cis-MR (drug-target): instruments restricted to +/- cis_window of each cis gene.
# This is the email's GOOD-example step (e.g. BCKDK/DBT/PPM1K -> BCAA -> T2D).
# snakemake@input$exposure / $outcome ; params: cis_window, r2, genes, panel
suppressMessages({library(data.table); library(MendelianRandomization)})
set.seed(snakemake@params$seed)   # reproducibility: global seed from config

exp <- fread(snakemake@input$exposure)
out <- fread(snakemake@input$outcome)
win <- snakemake@params$cis_window; r2 <- snakemake@params$r2
genes <- snakemake@params$genes;    panel <- snakemake@params$panel

# ---- TODO ----
# For each gene in `genes`:
#   1. Get gene coords (GRCh38, e.g. via biomaRt / a static bed).
#   2. Instruments = exposure SNPs within +/- win of the gene, pval < 5e-8, pruned r2 < r2.
#   3. Harmonize with outcome; IVW (weights="delta") + MR-Egger.
#   4. Row per gene: gene, nsnp, b, se, pval, egger_intercept_p.
# cis-MR is pleiotropy-robust (local instruments) but power-limited (few SNPs) -> MCQ.

res <- data.table(gene=unlist(genes), nsnp=NA_integer_, b=NA_real_, se=NA_real_, pval=NA_real_)
fwrite(res, snakemake@output[[1]], sep = "\t")
