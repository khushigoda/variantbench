#!/usr/bin/env Rscript
# Step 1 — fixed-effect IVW meta-analysis of the per-cohort normalized files.
# snakemake@input       : results/norm/<metab>.EstBB.norm.tsv.gz, .UKBB_EUR.norm.tsv.gz
# snakemake@output$meta : results/meta/<metab>.meta.tsv.gz
# snakemake@output$plot : results/meta/<metab>.concordance.png
suppressMessages({library(data.table); library(ggplot2)})

files <- snakemake@input
a <- fread(files[[1]]); b <- fread(files[[2]])

# ---- TODO ----
# 1. Merge a,b on variant_id (chr:pos:ref:alt already build-consistent from Step 0).
# 2. Allele-align: where ea/oa swapped, flip sign of beta and eaf; drop strand-ambiguous
#    A/T & G/C at MAF>0.4.
# 3. Fixed-effect IVW: w=1/se^2; beta_meta=sum(w*beta)/sum(w); se_meta=sqrt(1/sum(w)).
#    pval_meta from z=beta_meta/se_meta. Cochran's Q + het_p. direction string "+-".
# 4. Concordance plot: EstBB beta vs UKBB beta (expect tight diagonal).
# 5. VALIDATION (separate check, not this rule): correlate beta_meta vs published
#    meta_EUR accession -> r should be ~1. Note this in the report.

meta <- a  # <-- replace with merged/meta table incl beta_meta, se_meta, pval_meta, het_q, het_p, direction
fwrite(meta, snakemake@output$meta, sep = "\t", compress = "gzip")

p <- ggplot(meta, aes(beta, beta)) + geom_point(alpha=.3) + theme_bw() +
     labs(title = paste(snakemake@wildcards$metabolite, "cohort concordance"),
          x = "EstBB beta", y = "UKBB_EUR beta")
ggsave(snakemake@output$plot, p, width = 5, height = 5, dpi = 120)
