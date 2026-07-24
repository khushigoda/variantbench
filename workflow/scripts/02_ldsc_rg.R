#!/usr/bin/env Rscript
# Step 2 (rg) — run LDSC genetic-correlation across all metabolite pairs, parse to tsv.
# The h2 and munge steps are pure shell rules; this wrapper does --rg then tidies logs.
# snakemake@input  : results/ldsc/<metab>.sumstats.gz (all metabolites)
# snakemake@output : results/ldsc/rg.tsv
suppressMessages({library(data.table)})

sumstats <- unlist(snakemake@input)
ld       <- "data/ref/eur_w_ld_chr/"

# ---- TODO ----
# 1. Build the --rg comma list (first trait vs rest, or all pairwise).
# 2. system2("conda", c("run","-n","ldsc","python","~/ldsc/ldsc.py","--rg",
#      paste(sumstats, collapse=","), "--ref-ld-chr", ld, "--w-ld-chr", ld,
#      "--out","results/ldsc/rg")).
# 3. Parse the "Summary of Genetic Correlation Results" table from rg.log -> data.table.

rg <- data.table(p1=character(), p2=character(), rg=numeric(), se=numeric(), p=numeric())
fwrite(rg, snakemake@output[[1]], sep = "\t")
