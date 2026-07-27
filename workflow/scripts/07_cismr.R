#!/usr/bin/env Rscript
# =============================================================================
# 07_cismr.R — cis-Mendelian randomization (drug-target): metabolite -> disease.
# =============================================================================
# Faithful to the manuscript's cis-MR (Methods, "Cis-Mendelian randomization"):
# for each target gene, instruments = exposure variants in the +/-200 kb region
# around the GENE BODY with MAF > 1%, greedy LD-pruned to r2 < 0.01 (same greedy
# algorithm as genome-wide MR); estimator = multiplicative random-effects IVW-MR
# with weights = "delta", plus MR-Egger as the pleiotropy-robust companion.
# One causal estimate PER GENE.
#
# This differs from genome-wide MR (06_gwmr.R) in exactly three places:
#   (a) candidate instruments come from a per-gene cis-window slice of the
#       per-variant exposure meta stats (tabix), NOT the genome-wide lead table;
#   (b) IVW uses weights = "delta" (second-order SEs) — the paper specified this
#       for cis-MR (06 genome-wide MR used the default);
#   (c) MR-Egger is run per gene alongside IVW (its intercept flags directional
#       pleiotropy — the whole point of restricting to a cis region).
# The LD-matrix build, prune_by_LD(), and harmonization are reused verbatim.
#
# METHOD, per gene
#   1. Candidates = exposure meta variants within [start-win, end+win] with
#      MAF > MAF_MIN (eaf-derived) and P < GWS_P.
#   2. Greedy LD pruning: r2 from the 1000G-EUR panel via plink2 (one square r
#      matrix for the gene's chromosome slice), prune_by_LD() -> r2 < R2.
#   3. Harmonize the pruned instruments against {OUTCOME} on rsID (align to the
#      exposure effect allele, drop strand-ambiguous palindromes).
#   4. IVW-MR mr_ivw(model="random", weights="delta"); MR-Egger mr_egger().
#      Emit one row: gene, nsnp, IVW b/se/CI/P/OR, Egger b/P + intercept/P.
#
# INPUTS (env vars; self-bridged via `conda run -n variantbench-r` in the rule)
#   TRAIT      exposure/metabolite label            [LDL_C]
#   OUTCOME    outcome label (report/output naming)  [CAD]
#   GENES_TSV  cis-gene coord table (gene ensg chr start end, GRCh38 body)
#                                                    [config/cis_genes.{TRAIT}.tsv]
#   EXP_BGZ    exposure per-variant meta stats (tabix'd bgz)
#                                                    [results/meta/{TRAIT}.meta.tsv.bgz]
#   OUTCOME_GZ outcome GWAS sumstats (gz, harmonized) [data/raw/CAD.GCST90132314.h.tsv.gz]
#   REFDIR     1000G-EUR plink2 panel dir            [data/ref/1kg_phase3_hg38]
#   KEEP       EUR sample keep-file                  [{REFDIR}/1kg_eur.keep]
#   PLINK2     plink2 launcher                        [conda run -n plink2 plink2]
#   OUTDIR     output dir                            [results/mr]
#   CIS_WINDOW +/- bp around the gene body           [200000]
#   GWS_P      genome-wide significance              [5e-8]
#   MAF_MIN    instrument MAF floor                   [0.01]
#   R2         greedy-pruning LD r2 threshold         [0.01]
#   SEED       RNG seed (reproducibility)             [1123]
#
# OUTPUTS
#   {OUTDIR}/{TRAIT}.{OUTCOME}.cismr.tsv             one IVW+Egger row per gene
#   {OUTDIR}/{TRAIT}.{OUTCOME}.cismr.instruments.tsv per-SNP bx/by, all genes
#   {OUTDIR}/{TRAIT}.{OUTCOME}.cismr.png             per-gene IVW forest (OR + CI)
# =============================================================================
suppressMessages({ library(data.table); library(MendelianRandomization); library(ggplot2) })

env <- function(k, d="") Sys.getenv(k, unset=d)
TRAIT      <- env("TRAIT",   "LDL_C")
OUTCOME    <- env("OUTCOME", "CAD")
GENES_TSV  <- env("GENES_TSV",  sprintf("config/cis_genes.%s.tsv", TRAIT))
EXP_BGZ    <- env("EXP_BGZ",     sprintf("results/meta/%s.meta.tsv.bgz", TRAIT))
OUTCOME_GZ <- env("OUTCOME_GZ", "data/raw/CAD.GCST90132314.h.tsv.gz")
REFDIR     <- env("REFDIR",  "data/ref/1kg_phase3_hg38")
KEEP       <- env("KEEP",    file.path(REFDIR, "1kg_eur.keep"))
PLINK2     <- env("PLINK2",  "conda run -n plink2 plink2")
OUTDIR     <- env("OUTDIR",  "results/mr")
CIS_WINDOW <- as.integer(env("CIS_WINDOW", "200000"))
GWS_P      <- as.numeric(env("GWS_P",   "5e-8"))
MAF_MIN    <- as.numeric(env("MAF_MIN", "0.01"))
R2         <- as.numeric(env("R2",      "0.01"))
SEED       <- as.integer(env("SEED",    "1123"))
set.seed(SEED)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("== 07 cis-MR | %s -> %s | +/-%dkb  MAF>%.3g P<%.1e greedy-r2<%.3g  IVW(delta)+Egger ==\n",
            TRAIT, OUTCOME, CIS_WINDOW %/% 1000L, MAF_MIN, GWS_P, R2))

genes <- fread(GENES_TSV, skip = "gene\t")   # jump past leading '#' provenance comments
stopifnot(all(c("gene","chr","start","end") %in% names(genes)))
genes[, chr := as.character(chr)]
cat(sprintf("[0] cis genes: %s\n", paste(genes$gene, collapse = ", ")))

# ---- shared helper: exposure cis-window slice from the tabix'd meta stats -----
# meta header: SNP chr pos ea oa eaf beta se z pval neg_log10p n direction het_q het_p n_cohorts
exp_region <- function(chr, start, end) {
  reg <- sprintf("%s:%d-%d", sub("^chr","",as.character(chr)), start, end)
  txt <- tryCatch(system2("tabix", c(shQuote(EXP_BGZ), reg), stdout = TRUE, stderr = FALSE),
                  error = function(e) character(0))
  if (length(txt) == 0) return(data.table())
  d <- fread(text = paste(txt, collapse = "\n"), header = FALSE)
  setnames(d, c("SNP","chr","pos","ea","oa","eaf","beta","se","z","pval",
                "neg_log10p","n","direction","het_q","het_p","n_cohorts")[seq_len(ncol(d))])
  d[, maf := pmin(eaf, 1 - eaf)]
  d
}

# ---- shared helper: r2 among a chromosome's candidates from the panel ---------
# (identical build to 06_gwmr.R; here always a single-chromosome slice per gene)
ld_r2 <- function(ids, chr, tag) {
  n <- length(ids); r2mat <- matrix(0, n, n, dimnames = list(ids, ids)); diag(r2mat) <- 1
  if (n < 2) return(list(r2 = r2mat, matched = ids))
  pgen <- file.path(REFDIR, sprintf("chr%s_hg38.pgen",       chr))
  pvar <- file.path(REFDIR, sprintf("chr%s_hg38_rs.pvar.zst", chr))
  psam <- file.path(REFDIR, sprintf("chr%s_hg38.psam",       chr))
  if (!file.exists(pgen)) { cat(sprintf("   [warn] panel missing for chr%s — candidates stay independent\n", chr)); return(list(r2 = r2mat, matched = character(0))) }
  writeLines(ids, paste0(tag, ".ids"))
  cmd <- sprintf(paste("%s --pgen %s --pvar %s --psam %s --keep %s",
                       "--extract %s --r-unphased square ref-based --out %s"),
                 PLINK2, shQuote(pgen), shQuote(pvar), shQuote(psam), shQuote(KEEP),
                 shQuote(paste0(tag, ".ids")), shQuote(tag))
  st <- suppressWarnings(system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE))
  vcor <- paste0(tag, ".unphased.vcor1"); vars <- paste0(tag, ".unphased.vcor1.vars")
  if (st != 0 || !file.exists(vcor) || !file.exists(vars)) {
    cat(sprintf("   [warn] plink2 r failed for chr%s — candidates stay independent\n", chr))
    file.remove(Sys.glob(paste0(tag, "*"))); return(list(r2 = r2mat, matched = character(0)))
  }
  v  <- readLines(vars); rr <- as.matrix(fread(vcor, header = FALSE)); dimnames(rr) <- list(v, v)
  rr[is.na(rr)] <- 0
  common <- intersect(v, ids); r2mat[common, common] <- rr[common, common]^2
  # A variant monomorphic in the EUR panel gets r=NaN->0 above, zeroing its OWN
  # diagonal; prune_by_LD then never drops it (0 > r2 is FALSE) -> infinite loop.
  r2mat[!is.finite(r2mat)] <- 0; diag(r2mat) <- 1
  file.remove(Sys.glob(paste0(tag, "*")))
  list(r2 = r2mat, matched = common)
}

# ---- shared helper: greedy LD pruning (verbatim from authors' prune_by_LD) ----
prune_by_LD <- function(cand, r2mat, r2_limit) {
  s <- cand[order(cand$pval)]; keep <- character(0)
  while (nrow(s) > 0) {
    lead <- s$SNP[1]; keep <- c(keep, lead)
    in_ld <- names(which(r2mat[lead, ] > r2_limit))
    s <- s[!(SNP %in% in_ld)]
  }
  cand[SNP %in% keep]
}

# ---- shared helper: harmonize instruments against the outcome GWAS ------------
# Outcome GWAS is tabix-indexed on chr/pos; fetch ONLY each gene's pruned cis
# instruments (few-KB indexed read) instead of reading all ~19.5M rows upfront.
# Same build (GRCh38) both sides -> identical rows; rsID (col 26) kept for the
# rsID merge, so results match a full-file read exactly. If the .tbi is missing
# the pipeline builds it: `tabix -s 1 -b 2 -e 2 -S 1 <outcome>.h.tsv.gz`.
read_outcome_at <- function(pos_dt) {
  regs <- sprintf("%s:%d-%d", sub("^chr","",as.character(pos_dt$chr)), pos_dt$pos, pos_dt$pos)
  txt  <- tryCatch(system2("tabix", c(shQuote(OUTCOME_GZ), regs), stdout = TRUE, stderr = FALSE),
                   error = function(e) character(0))
  if (length(txt) == 0) return(data.table())
  d <- fread(text = paste(txt, collapse = "\n"), header = FALSE, sep = "\t")
  d <- d[, .(SNP = sub("\r$","", as.character(V26)),
             o_ea = V3, o_oa = V4, o_beta = as.numeric(V5),
             o_se = as.numeric(V6), o_pval = as.numeric(V8))]
  unique(d[SNP != "" & SNP != "NA"], by = "SNP")
}
harmonize <- function(instr) {
  oc <- read_outcome_at(instr)
  if (nrow(oc) == 0) return(data.table())
  h <- merge(instr[, .(SNP, ea, oa, e_beta = beta, e_se = se)], oc, by = "SNP")
  if (nrow(h) == 0) return(h)
  up <- function(x) toupper(as.character(x))
  h[, `:=`(ea = up(ea), oa = up(oa), o_ea = up(o_ea), o_oa = up(o_oa))]
  comp <- c(A="T", T="A", C="G", G="C")
  h <- h[!(o_ea == comp[ea] & o_oa == comp[oa] | ea == oa)]
  h <- h[(ea == o_ea & oa == o_oa) | (ea == o_oa & oa == o_ea)]
  h[, flip := as.integer(ea == o_oa & oa == o_ea)]
  h[, by_aligned := ifelse(flip == 1L, -o_beta, o_beta)]
  h
}

# ---- per-gene cis-MR ----------------------------------------------------------
res_rows <- list(); instr_rows <- list()
for (i in seq_len(nrow(genes))) {
  g <- genes[i]; sym <- g$gene; chr <- g$chr
  w0 <- max(1L, g$start - CIS_WINDOW); w1 <- g$end + CIS_WINDOW
  cand <- exp_region(chr, w0, w1)
  cand <- cand[maf > MAF_MIN & pval < GWS_P]
  cat(sprintf("[%s] window chr%s:%d-%d  candidates(MAF>%.3g,P<%.1e)=%d\n",
              sym, chr, w0, w1, MAF_MIN, GWS_P, nrow(cand)))
  if (nrow(cand) < 1) { cat(sprintf("   [%s] no cis instruments — skipped\n", sym)); next }

  ld <- ld_r2(cand$SNP, chr, file.path(OUTDIR, sprintf(".cisld_%s_%s", TRAIT, sym)))
  instr <- prune_by_LD(cand, ld$r2, R2)
  cat(sprintf("   [%s] LD coverage %d/%d; pruned -> %d instruments\n",
              sym, length(ld$matched), nrow(cand), nrow(instr)))

  h <- harmonize(instr)
  cat(sprintf("   [%s] harmonized instruments: %d\n", sym, nrow(h)))
  if (nrow(h) < 1) { cat(sprintf("   [%s] no harmonized instruments — skipped\n", sym)); next }
  h[, gene := sym]; instr_rows[[sym]] <- h

  # cis-MR needs >= 3 SNPs for MR-Egger (intercept + slope + residual df).
  mri <- mr_input(exposure = TRAIT, outcome = sym, snps = h$SNP,
                  bx = h$e_beta, bxse = h$e_se, by = h$by_aligned, byse = h$o_se)
  ivw <- mr_ivw(mri, model = "random", weights = "delta")     # paper: weights = "delta"
  eg  <- tryCatch(mr_egger(mri), error = function(e) NULL)     # >= 3 SNPs required
  res_rows[[sym]] <- data.table(
    gene = sym, chr = chr, nsnp = ivw@SNPs,
    ivw_b = ivw@Estimate, ivw_se = ivw@StdError,
    ivw_lo = ivw@CILower, ivw_up = ivw@CIUpper, ivw_p = ivw@Pvalue,
    ivw_or = exp(ivw@Estimate), ivw_or_lo = exp(ivw@CILower), ivw_or_up = exp(ivw@CIUpper),
    egger_b        = if (!is.null(eg)) eg@Estimate         else NA_real_,
    egger_p        = if (!is.null(eg)) eg@Pvalue.Est       else NA_real_,
    egger_intercept= if (!is.null(eg)) eg@Intercept        else NA_real_,
    egger_int_p    = if (!is.null(eg)) eg@Pvalue.Int       else NA_real_)
  cat(sprintf("   [%s] IVW OR=%.3f (%.3f,%.3f) P=%.2e | Egger int P=%s  nsnp=%d\n",
              sym, exp(ivw@Estimate), exp(ivw@CILower), exp(ivw@CIUpper), ivw@Pvalue,
              if (!is.null(eg)) sprintf("%.2g", eg@Pvalue.Int) else "NA(<3 SNP)", ivw@SNPs))
}
if (length(res_rows) == 0) stop("no gene yielded a cis-MR estimate")
res <- rbindlist(res_rows)
fwrite(res, file.path(OUTDIR, sprintf("%s.%s.cismr.tsv", TRAIT, OUTCOME)), sep = "\t")
fwrite(rbindlist(instr_rows)[, .(gene, SNP, ea, oa, bx = e_beta, bxse = e_se,
                                 by = by_aligned, byse = o_se, flip, o_pval)],
       file.path(OUTDIR, sprintf("%s.%s.cismr.instruments.tsv", TRAIT, OUTCOME)), sep = "\t")
cat(sprintf("[out] %d genes with cis-MR estimates\n", nrow(res)))

# ---- forest plot: per-gene IVW OR with 95% CI --------------------------------
res[, gene := factor(gene, levels = rev(gene))]
p <- ggplot(res, aes(ivw_or, gene)) +
  geom_vline(xintercept = 1, colour = "grey70", linewidth = 0.3, linetype = 2) +
  geom_errorbarh(aes(xmin = ivw_or_lo, xmax = ivw_or_up), height = 0.18, colour = "#2c3e50") +
  geom_point(size = 2.4, colour = "#c0392b") +
  labs(x = sprintf("%s -> %s  cis-MR OR per SD (95%% CI)", TRAIT, OUTCOME), y = NULL,
       title = sprintf("cis-MR (drug-target): %s -> %s", TRAIT, OUTCOME),
       subtitle = "multiplicative random-effects IVW (weights = delta)") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 10))
ggsave(file.path(OUTDIR, sprintf("%s.%s.cismr.png", TRAIT, OUTCOME)),
       plot = p, width = 5.4, height = 3.4, dpi = 150, bg = "white")
cat("[done] 07 cis-MR\n")
