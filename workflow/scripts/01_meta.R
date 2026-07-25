#!/usr/bin/env Rscript
# Step 1 — fixed-effect IVW meta-analysis, EstBB + UKBB_EUR -> meta.
#
# Reads the RAW GWAS-SSF files directly (column-selective fread; NO Step-0 normalized
# intermediate is written). NO-FLIP join: variants are matched on chr:pos:effect:other,
# so only identically-oriented alleles combine and single-cohort variants are retained
# (NA weight -> 0). This mirrors the authors' metaanalysis.R and lets us OBSERVE the cost
# of not flipping (reported as a concordance diagnostic) rather than pre-correcting.
#
#   input$est  : raw EstBB    .tsv.gz     input$ukbb : raw UKBB_EUR .tsv.gz
#   input$val  : raw published meta_EUR .tsv.gz  (validation TARGET, not a meta input)
#   output$meta    : results/meta/<metab>.meta.tsv.gz
#   output$plot    : results/meta/<metab>.concordance.png   (EstBB beta vs UKBB beta)
#   output$valplot : results/meta/<metab>.validation.png    (our beta vs published beta)
#   output$valtsv  : results/meta/<metab>.validation.tsv    (r, slope, n, join counts)
suppressMessages({library(data.table); library(ggplot2)})
set.seed(snakemake@params$seed)   # reproducibility: global seed from config
setDTthreads(1)                   # bound memory/CPU on the 8 GB box

# GWAS-SSF columns 01 actually consumes (confirm names against 00_inspect findings).
SEL <- c("chromosome", "base_pair_location", "effect_allele", "other_allele",
         "beta", "standard_error", "effect_allele_frequency", "rsid", "n")

# MEMORY NOTE (8 GB box): the join is on (chr,pos,ea,oa) as SEPARATE columns, NOT a
# pasted "chr:pos:ea:oa" string. A 26M-row character key vector is ~2 GB and was
# materialized in every table (est/ukbb/m/out/val) — the main cause of swap-thrash.
# Joining on the 4 columns is identical in result, adds zero new columns, and keeps
# chr/pos/ea/oa populated for every outer-join row (so no key-splitting recovery).
# chr is stored as INTEGER (1..22, X->23); non-autosomal/non-numeric contigs drop out
# (LDSC + this analysis are autosomal EUR anyway).
JKEY <- c("chr", "pos", "ea", "oa")

read_cohort <- function(f, tag, keep_rsid = FALSE) {
  message(sprintf("  [%s] reading %s ...", tag, basename(f)))
  d <- fread(f, select = SEL, showProgress = FALSE)
  setnames(d,
    c("chromosome", "base_pair_location", "effect_allele", "other_allele",
      "standard_error", "effect_allele_frequency"),
    c("chr", "pos", "ea", "oa", "se", "eaf"))
  d[, chr := toupper(sub("^chr", "", as.character(chr)))]
  d[chr == "X", chr := "23"]
  d[, chr := suppressWarnings(as.integer(chr))]
  d <- d[is.finite(beta) & is.finite(se) & se > 0 & !is.na(chr)]
  d <- unique(d, by = JKEY)
  # keep only what the downstream steps need, per cohort (drop early to save memory)
  if (keep_rsid) d <- d[, .(chr, pos, ea, oa, eaf, beta, se, n, rsid)]  # est: +eaf +rsid
  else           d <- d[, .(chr, pos, ea, oa, beta, se, n)]             # ukbb: minimal
  message(sprintf("  [%s] %d variants after QC", tag, nrow(d)))
  d
}

est  <- read_cohort(snakemake@input$est,  "EstBB",   keep_rsid = TRUE)
ukbb <- read_cohort(snakemake@input$ukbb, "UKBB_EUR")

# ---- no-flip full (outer) join on (chr,pos,ea,oa) ----
m <- merge(
  est[,  .(chr, pos, ea, oa, eaf, beta_e = beta, se_e = se, n_e = n)],
  ukbb[, .(chr, pos, ea, oa, beta_u = beta, se_u = se, n_u = n)],
  by = JKEY, all = TRUE)
rm(ukbb); invisible(gc())                     # ukbb done after the merge
est <- est[, .(chr, pos, ea, oa, rsid)]        # slim est to just the rsID lookup
invisible(gc())
message(sprintf("  merged: %d union variants", nrow(m)))

# per-cohort inverse-variance weights; absent cohort -> weight 0 (variant kept, single-cohort)
m[, `:=`(w_e = fifelse(is.finite(se_e), 1 / se_e^2, 0),
         w_u = fifelse(is.finite(se_u), 1 / se_u^2, 0),
         b_e = fifelse(is.finite(beta_e), beta_e, 0),
         b_u = fifelse(is.finite(beta_u), beta_u, 0))]
m[, w_sum := w_e + w_u]
m <- m[w_sum > 0]

# fixed-effect IVW
m[, `:=`(beta = (w_e * b_e + w_u * b_u) / w_sum,
         se   = sqrt(1 / w_sum),
         n    = fifelse(is.na(n_e), 0, as.numeric(n_e)) + fifelse(is.na(n_u), 0, as.numeric(n_u)),
         n_cohorts = (w_e > 0) + (w_u > 0))]
m[, z := beta / se]                                     # SIGNED z (sign carries direction)
# two-sided p and -log10 p from z; -log10 p in log-space to stay underflow-safe at huge z
m[, pval       := 2 * pnorm(-abs(z))]
m[, neg_log10p := -(pnorm(-abs(z), log.p = TRUE) + log(2)) / log(10)]
# 2-study Cochran's Q (1 df) — heterogeneity diagnostic only
m[, het_q := fifelse(n_cohorts == 2, w_e * (b_e - beta)^2 + w_u * (b_u - beta)^2, NA_real_)]
m[, het_p := fifelse(is.finite(het_q), pchisq(het_q, df = 1, lower.tail = FALSE), NA_real_)]
m[, direction := paste0(fifelse(w_e > 0, fifelse(b_e >= 0, "+", "-"), "?"),
                        fifelse(w_u > 0, fifelse(b_u >= 0, "+", "-"), "?"))]

# concordance-diagnostic counts + the both-cohort beta pairs (extract BEFORE freeing m)
diag <- data.table(
  n_total      = nrow(m),
  n_both       = nrow(m[n_cohorts == 2]),
  n_estbb_only = nrow(m[w_e > 0 & w_u == 0]),
  n_ukbb_only  = nrow(m[w_u > 0 & w_e == 0]))
both <- m[n_cohorts == 2, .(beta_e, beta_u)]

out <- m[, .(chr, pos, ea, oa, eaf, beta, se, z, pval, neg_log10p, n,
             direction, het_q, het_p, n_cohorts)]
rm(m); invisible(gc())                          # m no longer needed
# SNP = native rsID from the raw files (col `rsid`, ~98.7% rs#####) where the variant
# exists in EstBB, else a positional id as fallback (built only for the ~1.3% unmatched).
# LDSC munge merges to w_hm3.snplist BY rsID, so carrying rsid here lets Step 2 match HapMap3.
out[, SNP := NA_character_]
out[est, SNP := i.rsid, on = JKEY]
rm(est); invisible(gc())
out[is.na(SNP), SNP := paste(chr, pos, ea, oa, sep = ":")]
setcolorder(out, c("SNP", "chr", "pos", "ea", "oa", "eaf", "beta", "se", "z",
                   "pval", "neg_log10p", "n", "direction", "het_q", "het_p", "n_cohorts"))
setorder(out, chr, pos)
fwrite(out, snakemake@output$meta, sep = "\t", compress = "gzip")
message(sprintf("  wrote meta table: %d variants", nrow(out)))

# ---- concordance diagnostic: the cost of NOT flipping ----
p1 <- ggplot(both, aes(beta_e, beta_u)) +
  geom_point(alpha = .2, size = .5) +
  geom_abline(slope = 1, intercept = 0, colour = "red") +
  theme_bw() +
  labs(title = paste(snakemake@wildcards$metabolite, "cohort concordance"),
       subtitle = sprintf("shared (both-cohort) %s / %s total variants",
                          format(diag$n_both, big.mark = ","),
                          format(diag$n_total, big.mark = ",")),
       x = "EstBB beta", y = "UKBB_EUR beta")
ggsave(snakemake@output$plot, p1, width = 5, height = 5, dpi = 120)
rm(both); invisible(gc())

# ---- validation vs published meta_EUR (matched-orientation subset) ----
message("  reading published meta_EUR for validation ...")
val <- fread(snakemake@input$val,
             select = c("chromosome", "base_pair_location", "effect_allele",
                        "other_allele", "beta"), showProgress = FALSE)
setnames(val, c("chr", "pos", "ea", "oa", "beta_pub"))
val[, chr := toupper(sub("^chr", "", as.character(chr)))][chr == "X", chr := "23"]
val[, chr := suppressWarnings(as.integer(chr))]
val <- val[!is.na(chr)]
val <- unique(val, by = JKEY)
cmp <- merge(out[, .(chr, pos, ea, oa, beta)],
             val[, .(chr, pos, ea, oa, beta_pub)], by = JKEY)
rm(val); invisible(gc())
r  <- if (nrow(cmp) > 2) cor(cmp$beta, cmp$beta_pub) else NA_real_
sl <- if (nrow(cmp) > 2) unname(coef(lm(beta_pub ~ beta, cmp))[2]) else NA_real_
fwrite(data.table(metab = snakemake@wildcards$metabolite,
                  n_compared = nrow(cmp), pearson_r = r, slope = sl,
                  n_both = diag$n_both, n_estbb_only = diag$n_estbb_only,
                  n_ukbb_only = diag$n_ukbb_only),
       snakemake@output$valtsv, sep = "\t")
p2 <- ggplot(cmp, aes(beta, beta_pub)) +
  geom_point(alpha = .2, size = .5) +
  geom_abline(slope = 1, intercept = 0, colour = "red") +
  theme_bw() +
  labs(title = paste(snakemake@wildcards$metabolite, "vs published meta_EUR"),
       subtitle = sprintf("Pearson r = %.4f  (n = %s matched)",
                          r, format(nrow(cmp), big.mark = ",")),
       x = "our beta_meta", y = "published meta_EUR beta")
ggsave(snakemake@output$valplot, p2, width = 5, height = 5, dpi = 120)
cat(sprintf("%s: meta %d variants; validation r=%.4f on %d shared\n",
            snakemake@wildcards$metabolite, nrow(out), r, nrow(cmp)))
