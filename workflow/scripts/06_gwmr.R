#!/usr/bin/env Rscript
# Step 6 — genome-wide Mendelian randomization: metabolite (exposure) -> disease (outcome).
# snakemake@input$exposure / $outcome ; snakemake@output
suppressMessages({library(data.table); library(MendelianRandomization)})

exp <- fread(snakemake@input$exposure)
out <- fread(snakemake@input$outcome)
gws <- snakemake@params$gws_p; r2 <- snakemake@params$r2
maf <- snakemake@params$maf;   panel <- snakemake@params$panel

# ---- TODO ----
# 1. Instruments = exposure lead variants (pval < gws), MAF > maf, greedy-pruned r2 < r2
#    against `panel`.
# 2. Harmonize exposure/outcome on variant_id: align ea/oa, flip beta sign as needed,
#    drop palindromic-ambiguous SNPs.
# 3. mr_input(bx, bxse, by, byse); mr_ivw(model="random") + mr_egger + mr_median.
# 4. Emit method, estimate (log-OR), se, pval, nsnp, egger_intercept + p (pleiotropy).
res <- data.table(method=c("IVW","Egger","WMedian"), b=NA_real_, se=NA_real_, pval=NA_real_, nsnp=NA_integer_)
fwrite(res, snakemake@output$tsv, sep = "\t")

# ---- PLOTS (§7 report) ----
# scatter: per-SNP outcome effect vs exposure effect with IVW/Egger/median slopes.
# loo:     leave-one-out IVW estimates to spot a single driving instrument.
# (The cross-pair FOREST plot is assembled in the report from all *.gwmr.tsv, not here.)
# Scaffold writes valid empty PNGs so the declared outputs exist during development.
suppressMessages(library(ggplot2))
ggsave(snakemake@output$scatter, ggplot() + theme_void(), width = 5, height = 5, dpi = 120)
ggsave(snakemake@output$loo,     ggplot() + theme_void(), width = 5, height = 5, dpi = 120)
