#!/usr/bin/env Rscript
# Step 2 (rg) — run LDSC genetic-correlation across all metabolite pairs, parse to tsv.
# The h2 and munge steps are pure shell rules; this wrapper does --rg then tidies logs.
# snakemake@input  : results/ldsc/<metab>.sumstats.gz (all metabolites)
# snakemake@output : results/ldsc/rg.tsv
suppressMessages({library(data.table)})
set.seed(snakemake@params$seed)   # reproducibility: global seed from config

sumstats <- unlist(snakemake@input)
ld       <- "data/ref/eur_w_ld_chr/"

# ---- TODO ----
# 1. Build the --rg comma list (first trait vs rest, or all pairwise).
# 2. system2("conda", c("run","-n","ldsc","python","~/ldsc/ldsc.py","--rg",
#      paste(sumstats, collapse=","), "--ref-ld-chr", ld, "--w-ld-chr", ld,
#      "--out","results/ldsc/rg")).
# 3. Parse the "Summary of Genetic Correlation Results" table from rg.log -> data.table.
# 4. Also parse each <metab>.h2.log for "Total Observed scale h2: X (SE)".

suppressMessages(library(ggplot2))
rg <- data.table(p1=character(), p2=character(), rg=numeric(), se=numeric(), p=numeric())
fwrite(rg, snakemake@output$rg, sep = "\t")

# ---- PLOTS (§3 report) ----
# h2 bar +/- SE across metabolites, and rg heatmap. Scaffold writes valid empty PNGs so
# the rule's declared outputs all exist during incremental development; replace with real
# plots once the logs are parsed above.
# TODO h2_plot: geom_col(aes(metab, h2)) + geom_errorbar(h2-se, h2+se)
# TODO rg_plot: geom_tile(aes(p1, p2, fill=rg)) + scale_fill_gradient2()
ggsave(snakemake@output$h2_plot, ggplot() + theme_void(), width = 4, height = 4, dpi = 120)
ggsave(snakemake@output$rg_plot, ggplot() + theme_void(), width = 4, height = 4, dpi = 120)
