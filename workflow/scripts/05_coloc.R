#!/usr/bin/env Rscript
# Step 5 — colocalization (coloc.abf) between the metabolite and a molecular QTL at a locus.
# snakemake@input$meta ; snakemake@output ; params: eqtl dataset, priors, locus
suppressMessages({library(data.table); library(coloc)})

meta   <- fread(snakemake@input$meta)
eqtl   <- snakemake@params$eqtl
priors <- snakemake@params$priors
locus  <- snakemake@params$locus   # list(chr=, gene=)

# ---- TODO ----
# 1. Subset meta to the coloc locus (chr == locus$chr, +/- window around the gene).
# 2. Pull matching eQTL Catalogue sumstats for `eqtl` dataset at the same locus
#    (FTP: ftp.ebi.ac.uk/pub/databases/spot/eQTL/), align variants.
# 3. Build coloc datasets (type "quant"): beta, varbeta=se^2, snp, MAF, N.
# 4. coloc.abf(d1, d2, p1=priors$p1, p2=priors$p2, p12=priors$p12); capture PP.H0..H4.
# 5. Report PP.H4 (shared causal variant). Mention CLPP-on-credible-sets as paper's method.

res <- data.table(metabolite=snakemake@wildcards$metabolite, gene=locus$gene,
                 PP.H0=NA_real_, PP.H1=NA_real_, PP.H2=NA_real_, PP.H3=NA_real_, PP.H4=NA_real_)
fwrite(res, snakemake@output[[1]], sep = "\t")
