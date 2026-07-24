#!/usr/bin/env Rscript
# Step 0 — normalize a raw sumstats file to the canonical "spine".
# Invoked by rules normalize_metabolite and normalize_outcome.
# snakemake@input[[1]]  : raw file. Builds by source:
#     metabolites .tsv.gz          -> GRCh38 (no liftover)
#     CAD .h.tsv.gz (harmonised)   -> GRCh38 (no liftover)
#     T2D .tsv.gz (DIAGRAM/hg19)    -> GRCh37 (Step 0 lifts to GRCh38)
# snakemake@output[[1]] : results/norm/<dataset>.norm.tsv.gz
# snakemake@params$build: source genome build
suppressMessages({library(data.table)})
set.seed(snakemake@params$seed)   # reproducibility: global seed from config

infile  <- snakemake@input[[1]]
outfile <- snakemake@output[[1]]
build   <- snakemake@params$build

dt <- fread(infile)

# ---- TODO: map source columns -> spine ----
# Spine columns: variant_id, chr, pos, ea, oa, eaf, beta, se, pval, n, info
# 1. Rename source cols (GWAS-SSF: chromosome, base_pair_location, effect_allele,
#    other_allele, effect_allele_frequency, beta, standard_error, ...).
# 2. Metabolites store -log10 p: pval <- 10^(-neg_log_10_p_value).
# 3. If build == "GRCh37": liftOver pos to GRCh38 (rtracklayer::liftOver + chain file),
#    drop unmapped, then continue.
# 4. Rebuild variant_id = paste(chr, pos, ref, alt) on GRCh38 (do NOT trust source rsID).
# 5. QC: metabolites INFO >= snakemake@params$info_min, MAC >= 20.
# 6. Keep spine columns in fixed order.

spine <- dt   # <-- replace with the transformed table
fwrite(spine, outfile, sep = "\t", compress = "gzip")
cat("normalized", infile, "->", outfile, "build", build, "\n")
