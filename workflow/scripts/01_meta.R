#!/usr/bin/env Rscript
# Step 1 — GENOME-WIDE fixed-effect IVW meta-analysis, EstBB + UKBB_EUR -> meta.
# Chromosome-streamed so it fits an 8 GB machine while keeping EVERY variant
# (no HapMap3 restriction) — the genome-wide table downstream steps
# (03 lead variants / Manhattan, 04 fine-mapping, 06/07 MR) all require.
#
# WHY STREAMING: a fixed-effect meta needs both cohorts' rows for a variant
# resident together. EstBB is 26M variants, UKBB_EUR is 96M; holding both full
# tables + their ~100M-row union does not fit in 8 GB (swap-thrash / OOM). Even
# the authors' metaanalysis.R accumulates the full genome-wide table before
# writing — fine on a big node, not here. So we mirror THEIR memory device
# (their Arrow open_dataset %>% filter(CHROM==c) reads one chromosome at a time)
# but go one step further: we also WRITE each chromosome's result incrementally
# (append to the gzip) and reduce validation to running SUFFICIENT STATISTICS,
# so peak memory is ONE chromosome, never the genome.
#
# NO-FLIP join (matches the authors): variants matched on chr:pos:effect:other,
# so only identically-oriented alleles combine; single-cohort variants retained
# (absent-cohort weight -> 0). We OBSERVE the cost of not flipping as a
# concordance diagnostic rather than pre-correcting.
#
#   input$est_dir  : dir of per-chr EstBB    split files  (chrN.tsv.gz)
#   input$ukbb_dir : dir of per-chr UKBB_EUR split files
#   input$val_dir  : dir of per-chr published meta_EUR split files (validation)
#   output$meta    : results/meta/<metab>.meta.tsv.gz      (GENOME-WIDE)
#   output$plot    : results/meta/<metab>.concordance.png  (EstBB beta vs UKBB beta)
#   output$valplot : results/meta/<metab>.validation.png   (meta_rep vs meta_og beta)
#   output$valtsv  : results/meta/<metab>.validation.tsv   (r, slope, n, counts)
suppressMessages({library(data.table); library(ggplot2)})
set.seed(snakemake@params$seed)   # reproducibility: global seed from config
setDTthreads(1)                   # bound memory/CPU on the 8 GB box

est_dir  <- snakemake@input$est_dir
ukbb_dir <- snakemake@input$ukbb_dir
val_dir  <- snakemake@input$val_dir
out_meta <- snakemake@output$meta
metab    <- snakemake@wildcards$metabolite

JKEY <- c("chr", "pos", "ea", "oa")

# ---- read one chromosome's split file, standardize names, QC ----
# Split files carry: chromosome base_pair_location effect_allele other_allele
#                    beta standard_error effect_allele_frequency [rsid] n
read_chr <- function(dir, chr_tok, chr_int) {
  f <- file.path(dir, paste0("chr", chr_tok, ".tsv.gz"))
  if (!file.exists(f)) return(NULL)
  d <- fread(f, showProgress = FALSE)
  ren <- c(chromosome = "chr", base_pair_location = "pos", effect_allele = "ea",
           other_allele = "oa", standard_error = "se", effect_allele_frequency = "eaf")
  setnames(d, names(ren), unname(ren), skip_absent = TRUE)
  # fread infers an all-NA column's type as LOGICAL. On chrX one cohort's rsid can be
  # entirely empty -> logical, which then breaks the character rsid coalesce. Force it
  # (and the allele cols, same failure mode) to the type the rest of the script expects.
  if ("rsid" %in% names(d)) d[, rsid := as.character(rsid)]
  if ("eaf"  %in% names(d)) d[, eaf  := as.numeric(eaf)]   # guard: all-NA chr -> logical via fread
  d[, `:=`(ea = as.character(ea), oa = as.character(oa))]
  d <- d[is.finite(beta) & is.finite(se) & se > 0]
  d[, chr := chr_int]                       # integer chr for this file
  d <- unique(d, by = c("pos", "ea", "oa"))
  d
}

# ---- discover chromosomes present (union of est + ukbb split dirs) ----
tok_of <- function(dir) sub("\\.tsv\\.gz$", "", sub("^chr", "",
              list.files(dir, pattern = "^chr.*\\.tsv\\.gz$")))
toks <- union(tok_of(est_dir), tok_of(ukbb_dir))
# map raw token -> integer chr (1..22, X->23); drop other contigs (autosomal+X EUR)
tok2int <- function(t) { u <- toupper(t); if (u == "X") 23L else suppressWarnings(as.integer(u)) }
toks <- toks[!is.na(vapply(toks, tok2int, integer(1)))]
toks <- toks[order(vapply(toks, tok2int, integer(1)))]
message(sprintf("  streaming %d chromosomes: %s", length(toks),
                paste(toks, collapse = ",")))

# ---- running accumulators (all O(1) memory) ----
vs <- c(n = 0, sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0)  # validation suff. stats
tot <- c(n_total = 0, n_both = 0, n_est = 0, n_ukbb = 0)
SPC <- 4000L                                # sampled points/chr for the scatter plots
conc_smp <- vector("list", length(toks))    # concordance: (beta_e, beta_u) shared
val_smp  <- vector("list", length(toks))    # validation:  (beta, beta_pub) matched
if (file.exists(out_meta)) file.remove(out_meta)   # fresh append target
first <- TRUE

for (ci in seq_along(toks)) {
  chr_tok <- toks[ci]; chr_int <- tok2int(chr_tok)
  est  <- read_chr(est_dir,  chr_tok, chr_int)
  ukbb <- read_chr(ukbb_dir, chr_tok, chr_int)
  if (is.null(est) && is.null(ukbb)) next
  if (is.null(est))  est  <- ukbb[0]
  if (is.null(ukbb)) ukbb <- est[0]

  # no-flip outer join on (chr,pos,ea,oa); rsid carried from both sides
  m <- merge(
    est[,  .(chr, pos, ea, oa, eaf_e = eaf, beta_e = beta, se_e = se, n_e = n, rsid_e = rsid)],
    ukbb[, .(chr, pos, ea, oa, eaf_u = eaf, beta_u = beta, se_u = se, n_u = n, rsid_u = rsid)],
    by = JKEY, all = TRUE)
  rm(est, ukbb)

  m[, `:=`(w_e = fifelse(is.finite(se_e), 1 / se_e^2, 0),
           w_u = fifelse(is.finite(se_u), 1 / se_u^2, 0),
           b_e = fifelse(is.finite(beta_e), beta_e, 0),
           b_u = fifelse(is.finite(beta_u), beta_u, 0))]
  m[, w_sum := w_e + w_u]
  m <- m[w_sum > 0]
  m[, `:=`(beta = (w_e * b_e + w_u * b_u) / w_sum,
           se   = sqrt(1 / w_sum),
           n    = fifelse(is.na(n_e), 0, as.numeric(n_e)) + fifelse(is.na(n_u), 0, as.numeric(n_u)),
           n_cohorts = (w_e > 0) + (w_u > 0))]
  # Combine effect-allele frequency across cohorts. The join is no-flip on (chr,pos,ea,oa),
  # so eaf_e and eaf_u are on the SAME effect allele and pool directly (like beta). Two-cohort
  # variants get the N-weighted mean; single-cohort variants keep their one cohort's eaf via
  # fcoalesce fall-through; a variant absent from both is genuinely NA. This populates eaf for
  # UKBB-only variants (previously NA because only EstBB's eaf was carried), which the Step-3
  # dual MAF threshold and Step-4 fine-mapping require.
  m[, eaf := {
       ne_ <- fifelse(is.na(n_e), 0, as.numeric(n_e))
       nu_ <- fifelse(is.na(n_u), 0, as.numeric(n_u))
       fcoalesce(
         fifelse(is.finite(eaf_e) & is.finite(eaf_u) & (ne_ + nu_) > 0,
                 (ne_ * eaf_e + nu_ * eaf_u) / (ne_ + nu_), NA_real_),
         eaf_e, eaf_u)
     }]
  m[, z := beta / se]                                   # SIGNED z (sign = direction)
  m[, pval       := 2 * pnorm(-abs(z))]
  m[, neg_log10p := -(pnorm(-abs(z), log.p = TRUE) + log(2)) / log(10)]
  m[, het_q := fifelse(n_cohorts == 2, w_e * (b_e - beta)^2 + w_u * (b_u - beta)^2, NA_real_)]
  m[, het_p := fifelse(is.finite(het_q), pchisq(het_q, df = 1, lower.tail = FALSE), NA_real_)]
  m[, direction := paste0(fifelse(w_e > 0, fifelse(b_e >= 0, "+", "-"), "?"),
                          fifelse(w_u > 0, fifelse(b_u >= 0, "+", "-"), "?"))]
  # coalesced native rsID (both cols forced to character in read_chr, so types match)
  m[, SNP := fcoalesce(as.character(rsid_e), as.character(rsid_u))]
  # both cohorts missing an rsID (e.g. some chrX variants) -> positional id, never NA
  m[is.na(SNP) | SNP == "", SNP := paste(chr, pos, ea, oa, sep = ":")]

  # running counts + concordance sample (both-cohort betas)
  tot["n_total"] <- tot["n_total"] + nrow(m)
  nb <- m[n_cohorts == 2, .N]; tot["n_both"] <- tot["n_both"] + nb
  tot["n_est"]  <- tot["n_est"]  + m[w_e > 0 & w_u == 0, .N]
  tot["n_ukbb"] <- tot["n_ukbb"] + m[w_u > 0 & w_e == 0, .N]
  both <- m[n_cohorts == 2, .(beta_e, beta_u)]
  if (nrow(both)) conc_smp[[ci]] <- both[sample(.N, min(.N, SPC))]

  out <- m[, .(SNP, chr, pos, ea, oa, eaf, beta, se, z, pval, neg_log10p, n,
               direction, het_q, het_p, n_cohorts)]
  setorder(out, pos)
  # INCREMENTAL WRITE: header on first chr only, append (concatenated gzip) after
  fwrite(out, out_meta, sep = "\t", compress = "gzip",
         append = !first, col.names = first)
  first <- FALSE

  # ---- validation for this chromosome: merge vs published, update suff. stats ----
  val <- read_chr(val_dir, chr_tok, chr_int)   # val has no rsid col
  if (!is.null(val)) {
    setnames(val, "beta", "beta_pub")
    cmp <- merge(out[, .(chr, pos, ea, oa, beta)],
                 val[, .(chr, pos, ea, oa, beta_pub)], by = JKEY)
    if (nrow(cmp)) {
      x <- cmp$beta; y <- cmp$beta_pub
      vs["n"]  <- vs["n"]  + length(x)
      vs["sx"] <- vs["sx"] + sum(x);  vs["sy"]  <- vs["sy"]  + sum(y)
      vs["sxx"]<- vs["sxx"]+ sum(x*x);vs["syy"] <- vs["syy"] + sum(y*y)
      vs["sxy"]<- vs["sxy"]+ sum(x*y)
      val_smp[[ci]] <- cmp[sample(.N, min(.N, SPC)), .(beta, beta_pub)]
    }
    rm(val, cmp)
  }
  rm(m, out, both); invisible(gc())
  message(sprintf("  chr%-3s done | cumulative meta rows: %s", chr_tok,
                  format(tot["n_total"], big.mark = ",")))
}
message(sprintf("  wrote GENOME-WIDE meta: %s variants -> %s",
                format(tot["n_total"], big.mark = ","), out_meta))

# ---- validation r + OLS slope from running sufficient statistics ----
n <- vs["n"]
r  <- if (n > 2) unname((n*vs["sxy"] - vs["sx"]*vs["sy"]) /
        sqrt((n*vs["sxx"] - vs["sx"]^2) * (n*vs["syy"] - vs["sy"]^2))) else NA_real_
sl <- if (n > 2) unname((n*vs["sxy"] - vs["sx"]*vs["sy"]) /
        (n*vs["sxx"] - vs["sx"]^2)) else NA_real_          # slope of beta_pub ~ beta
fwrite(data.table(metab = metab, n_compared = as.integer(n),
                  pearson_r = r, slope = sl, n_variants = as.integer(tot["n_total"]),
                  n_both = as.integer(tot["n_both"]),
                  n_estbb_only = as.integer(tot["n_est"]),
                  n_ukbb_only = as.integer(tot["n_ukbb"])),
       snakemake@output$valtsv, sep = "\t")

# ---- plots from the per-chr samples (a scatter of 100M points is unreadable) ----
conc <- rbindlist(conc_smp); vsc <- rbindlist(val_smp)
p1 <- ggplot(conc, aes(beta_e, beta_u)) +
  geom_point(alpha = .2, size = .5) +
  geom_abline(slope = 1, intercept = 0, colour = "red") + theme_bw() +
  labs(title = metab,
       x = "EstBB beta", y = "UKBB_EUR beta")
ggsave(snakemake@output$plot, p1, width = 5, height = 5, dpi = 120)
p2 <- ggplot(vsc, aes(beta, beta_pub)) +
  geom_point(alpha = .2, size = .5) +
  geom_abline(slope = 1, intercept = 0, colour = "red") + theme_bw() +
  labs(title = metab,
       x = "meta_rep beta", y = "meta_og beta")
ggsave(snakemake@output$valplot, p2, width = 5, height = 5, dpi = 120)
cat(sprintf("%s: GENOME-WIDE meta %s variants; validation r=%.4f slope=%.3f on %s shared\n",
            metab, format(tot["n_total"], big.mark = ","), r, sl,
            format(as.integer(n), big.mark = ",")))
