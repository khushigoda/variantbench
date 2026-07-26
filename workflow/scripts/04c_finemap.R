#!/usr/bin/env Rscript
# =============================================================================
# 04c_finemap.R — Step 4c: SuSiE fine-mapping (z-score mode), one locus at a time
# =============================================================================
# Consumes ONLY the per-locus Step-4 (make_ld) outputs — never re-slices the meta:
#   results/ld/<TRAIT>.loci.tsv                 (loop driver: gene chr lead_snp lead_pos ea oa z nlp start end)
#   results/ld/<TRAIT>.<gene>.ld                (signed-r matrix, tab, no header, N x N)
#   results/ld/<TRAIT>.<gene>.ld.vars           (variant IDs, 1/line = row/col order of .ld)
#   results/ld/<TRAIT>.<gene>.z.tsv             (ID POS z, SAME order as .ld.vars — make_ld's alignz guarantee)
#   results/loci/<TRAIT>.lead.tsv               (source of representative n per locus)
#
# Writes (single source of truth for all downstream steps — 04d annotation, plots, QC, coloc):
#   results/finemap/<TRAIT>.finemap.tsv         (per-variant: PIP, CS membership, r2-to-lead)
#   results/finemap/<TRAIT>.cs_summary.tsv      (per-credible-set headline table)
#   results/finemap/<TRAIT>.susie_qc.tsv        (per-locus convergence QC; kriging_s deferred)
#   results/finemap/<TRAIT>.<gene>.susie.rds    (full fit object — downstream reloads, NEVER re-fits)
#
# Method: susie_rss z-mode. Params mirror the authors (Tambets et al.) where they apply:
#   L=10, estimate_residual_variance=FALSE (out-of-sample LD, susieR issue #162),
#   estimate_prior_variance=TRUE, scaled_prior_variance=0.1, CS purity min_abs_corr=0.5.
# n = the lead variant's GWAS N (per-variant N is ~flat across a 1 Mb window, so lead N ~= max N;
#     matches the authors' max(N) choice while staying self-contained — no meta re-read).
# DEVIATION (documented): LD is OUT-OF-SAMPLE 1000G-EUR (525 founders), not in-sample UKBB.
#     Credible sets can be miscalibrated; the kriging_rss/estimate_s_rss LD-mismatch diagnostic
#     is deferred to a later QC step (kriging_s column is a stub here).
#
# Config is entirely env-var driven (nothing trait-hardcoded) so Snakemake wildcards and future
# traits/loci need only a config edit. Optional LOCI=gene[,gene] subsets for smoke tests.
# =============================================================================

suppressMessages({library(data.table); library(susieR)})

## ---- configuration (env-driven) --------------------------------------------
TRAIT   <- Sys.getenv("TRAIT",  "LDL_C")
SEED    <- as.integer(Sys.getenv("SEED", "1123"))
LDDIR   <- Sys.getenv("LDDIR",  "results/ld")
OUTDIR  <- Sys.getenv("OUTDIR", "results/finemap")
LEADTSV <- Sys.getenv("LEADTSV", file.path("results/loci", paste0(TRAIT, ".lead.tsv")))
LOCI    <- Sys.getenv("LOCI",   "")            # optional CSV gene subset (smoke test)
L_MAX   <- as.integer(Sys.getenv("L_MAX", "10"))
COV1    <- 0.95                                 # primary credible-set coverage
COV2    <- 0.99                                 # secondary coverage (also emitted)
MINCORR <- 0.5                                  # susie_get_cs purity floor (authors' min_cs_corr)
LOWPUR  <- 0.5                                  # low-purity flag threshold
GWS_P   <- as.numeric(Sys.getenv("GWS_P", "5e-8"))       # genome-wide sig (config params:gws_p)
ZTHRESH <- qnorm(GWS_P / 2, lower.tail = FALSE)          # |z| cutoff for a "real" CS (~5.45 at 5e-8)

set.seed(SEED)
suppressWarnings(setDTthreads(0))               # all cores for fread / matrix coercion
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("== 04c SuSiE fine-mapping | trait=%s seed=%d L=%d coverage=%.2f/%.2f ==\n",
            TRAIT, SEED, L_MAX, COV1, COV2))

## ---- loop driver + representative n ----------------------------------------
loci_path <- file.path(LDDIR, paste0(TRAIT, ".loci.tsv"))
stopifnot(file.exists(loci_path))
loci <- fread(loci_path)
if (nzchar(LOCI)) {
  want <- trimws(strsplit(LOCI, ",")[[1]])
  loci <- loci[gene %in% want]
  cat(sprintf("   LOCI subset -> %s (%d of requested %d found)\n",
              paste(loci$gene, collapse=","), nrow(loci), length(want)))
}
stopifnot(nrow(loci) > 0)

lead <- fread(LEADTSV)                            # has SNP + n columns
stopifnot(all(c("SNP","n") %in% names(lead)))
setkey(lead, SNP)

## ---- per-variant CS membership helper --------------------------------------
cs_membership <- function(cs_obj, p) {
  out <- rep(0L, p)
  if (length(cs_obj$cs)) for (k in seq_along(cs_obj$cs)) out[cs_obj$cs[[k]]] <- cs_obj$cs_index[k]
  out
}

## ---- fit one locus ---------------------------------------------------------
fit_one <- function(g) {
  row  <- loci[gene == g]
  lsnp <- row$lead_snp; lpos <- row$lead_pos
  base <- file.path(LDDIR, paste0(TRAIT, ".", g))

  vars <- fread(paste0(base, ".ld.vars"), header = FALSE)$V1
  zdt  <- fread(paste0(base, ".z.tsv"))                        # ID POS z
  R    <- as.matrix(fread(paste0(base, ".ld"), header = FALSE))
  R    <- (R + t(R)) / 2          # enforce exact symmetry (no-op numerically; silences susieR's internal warning)

  # make_ld guarantees identical ordering across the three files — assert it, don't assume it
  p <- length(vars)
  stopifnot(nrow(R) == p, ncol(R) == p, nrow(zdt) == p, all(zdt$ID == vars))

  # representative n = lead variant's GWAS N (fallback: max lead N if lead absent from table)
  n_loc <- lead[.(lsnp)]$n
  if (length(n_loc) == 0 || is.na(n_loc[1])) n_loc <- max(lead$n, na.rm = TRUE)
  n_loc <- as.numeric(n_loc[1])

  z <- zdt$z
  t0 <- Sys.time()
  fit <- susie_rss(z = z, R = R, n = n_loc, L = L_MAX,
                   estimate_residual_variance = FALSE,
                   estimate_prior_variance    = TRUE,
                   scaled_prior_variance       = 0.1)
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  cs95 <- susie_get_cs(fit, coverage = COV1, Xcorr = R, min_abs_corr = MINCORR)
  cs99 <- susie_get_cs(fit, coverage = COV2, Xcorr = R, min_abs_corr = MINCORR)
  pip  <- fit$pip

  # LD-mismatch diagnostic (kriging): s in [0,1], 0=consistent, ->1=severe out-of-sample mismatch
  s_est <- tryCatch(estimate_s_rss(z = z, R = R, n = n_loc),
                    error = function(e) NA_real_)

  lead_idx <- match(lsnp, vars)
  r2_lead  <- if (!is.na(lead_idx)) as.numeric(R[lead_idx, ])^2 else rep(NA_real_, p)

  finemap <- data.table(
    trait = TRAIT, gene = g, chr = row$chr,
    ID = vars, POS = zdt$POS, z = z,
    pip = pip, r2_lead = r2_lead,
    cs95 = cs_membership(cs95, p), cs99 = cs_membership(cs99, p),
    is_lead = vars == lsnp
  )

  cs_summary <- NULL
  if (length(cs95$cs)) {
    cs_summary <- rbindlist(lapply(seq_along(cs95$cs), function(k) {
      idx <- cs95$cs[[k]]; top <- idx[which.max(pip[idx])]
      cs_max_absz <- max(abs(z[idx]))
      data.table(
        trait = TRAIT, gene = g, chr = row$chr,
        cs = cs95$cs_index[k], cs_size = length(idx),
        top_variant = vars[top], top_pos = zdt$POS[top], top_pip = pip[top],
        cs_pip_sum = sum(pip[idx]),
        min_r2  = cs95$purity$min.abs.corr[k]^2,
        mean_r2 = cs95$purity$mean.abs.corr[k]^2,
        low_purity = cs95$purity$min.abs.corr[k] < LOWPUR,
        cs_max_absz = cs_max_absz,
        pass_filter = cs_max_absz > ZTHRESH,   # TRUE = contains a genome-wide-sig variant (real signal)
        coverage = COV1
      )
    }))
  }

  qc <- data.table(
    trait = TRAIT, gene = g, chr = row$chr,
    lead_snp = lsnp, lead_pos = lpos,
    n_variants = p, n_used = n_loc,
    converged = fit$converged, niter = fit$niter,
    n_cs95 = length(cs95$cs), n_cs99 = length(cs99$cs),
    n_cs95_pass = if (length(cs95$cs)) sum(sapply(cs95$cs, function(i) max(abs(z[i])) > ZTHRESH)) else 0L,
    max_pip = max(pip), elbo = tail(fit$elbo, 1),
    kriging_s = s_est,                             # LD-mismatch: 0=consistent, ->1=severe out-of-sample mismatch
    elapsed_s = round(el, 2)
  )

  saveRDS(fit, file.path(OUTDIR, paste0(TRAIT, ".", g, ".susie.rds")))
  list(finemap = finemap, cs = cs_summary, qc = qc)
}

## ---- loop (robust: one bad locus cannot kill the run) ----------------------
all_fm <- list(); all_cs <- list(); all_qc <- list()
for (g in loci$gene) {
  cat(sprintf("  -- %s  (lead %s) --\n", g, loci[gene == g]$lead_snp))
  res <- tryCatch(fit_one(g), error = function(e) { cat("     [FAIL] ", conditionMessage(e), "\n"); NULL })
  if (is.null(res)) {
    r <- loci[gene == g]
    all_qc[[g]] <- data.table(trait=TRAIT, gene=g, chr=r$chr, lead_snp=r$lead_snp, lead_pos=r$lead_pos,
                              n_variants=NA_integer_, n_used=NA_real_, converged=FALSE, niter=NA_integer_,
                              n_cs95=NA_integer_, n_cs99=NA_integer_, max_pip=NA_real_, elbo=NA_real_,
                              kriging_s=NA_real_, elapsed_s=NA_real_)
    next
  }
  all_fm[[g]] <- res$finemap; all_cs[[g]] <- res$cs; all_qc[[g]] <- res$qc
  cat(sprintf("     converged=%s niter=%s n_cs95=%d n_cs99=%d max_pip=%.3f  (%.1fs)\n",
              res$qc$converged, res$qc$niter, res$qc$n_cs95, res$qc$n_cs99, res$qc$max_pip, res$qc$elapsed_s))
}

## ---- write outputs (header-only guards keep the schema stable if a table is empty) ----
fm_out <- rbindlist(all_fm, fill = TRUE)
cs_out <- rbindlist(all_cs, fill = TRUE)
qc_out <- rbindlist(all_qc, fill = TRUE)

fwrite(fm_out, file.path(OUTDIR, paste0(TRAIT, ".finemap.tsv")),    sep = "\t")
if (nrow(cs_out) == 0)
  cs_out <- data.table(trait=character(), gene=character(), chr=integer(), cs=integer(), cs_size=integer(),
                       top_variant=character(), top_pos=integer(), top_pip=double(), cs_pip_sum=double(),
                       min_r2=double(), mean_r2=double(), low_purity=logical(),
                       cs_max_absz=double(), pass_filter=logical(), coverage=double())
fwrite(cs_out, file.path(OUTDIR, paste0(TRAIT, ".cs_summary.tsv")), sep = "\t")
# filtered view: only CS containing a genome-wide-significant variant (drops LD-mismatch artifacts)
cs_filt <- cs_out[pass_filter == TRUE]
fwrite(cs_filt, file.path(OUTDIR, paste0(TRAIT, ".cs_summary.filtered.tsv")), sep = "\t")
fwrite(qc_out, file.path(OUTDIR, paste0(TRAIT, ".susie_qc.tsv")),   sep = "\t")

cat(sprintf("\n[done] %d loci fit | %d CS (95%%, all) | %d pass filter | %d variants total\n",
            nrow(qc_out), nrow(cs_out), nrow(cs_filt), nrow(fm_out)))
cat(sprintf("       -> %s/{%s.finemap.tsv, %s.cs_summary.tsv, %s.cs_summary.filtered.tsv, %s.susie_qc.tsv} + per-locus .susie.rds\n",
            OUTDIR, TRAIT, TRAIT, TRAIT, TRAIT))
