#!/usr/bin/env Rscript
# Step 3 — genome-wide-significant lead variants + locus definition, Manhattan plot.
# snakemake@input$meta   : results/meta/<metab>.meta.tsv.gz
# snakemake@output$lead  : results/loci/<metab>.lead.tsv
# snakemake@output$plot  : results/loci/<metab>.manhattan.png
suppressMessages({library(data.table); library(ggplot2)})

meta  <- fread(snakemake@input$meta)
gws   <- snakemake@params$gws_p
win   <- snakemake@params$window
r2    <- snakemake@params$r2
panel <- snakemake@params$panel

# ---- TODO ----
# 1. Filter pval_meta < gws.
# 2. LD-aware clumping: PLINK --clump on `panel` (r2 threshold), or manual: sort by p,
#    greedily take top SNP, remove all within +/- win/2 OR r2 >= clump_r2_lead, repeat.
# 3. Emit lead.tsv: variant_id, chr, pos, ea, oa, beta, se, pval, locus_start, locus_end.
# 4. Manhattan coloured by chr with lead loci labelled.

lead <- meta[pval < gws][order(pval)]
fwrite(lead, snakemake@output$lead, sep = "\t")

p <- ggplot(meta, aes(pos, -log10(pval))) + geom_point(size=.3) +
     geom_hline(yintercept=-log10(gws), linetype=2, colour="red") + theme_bw() +
     labs(title = paste(snakemake@wildcards$metabolite, "Manhattan"))
ggsave(snakemake@output$plot, p, width = 9, height = 4, dpi = 120)
