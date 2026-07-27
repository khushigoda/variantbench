#!/usr/bin/env Rscript
# =============================================================================
# 06_gwmr.R — genome-wide Mendelian randomization: metabolite -> disease.
# =============================================================================
# Faithful to the manuscript's genome-wide MR (Methods, "Genome-wide Mendelian
# randomization"): instruments = the trait's lead variants with MAF > 1%, greedy
# LD-pruned to r2 < 0.01; estimator = multiplicative random-effects IVW-MR from
# the MendelianRandomization R package, as recommended by Burgess et al. 2023.
# A SINGLE estimator, one causal estimate. (Robust/sensitivity methods — Egger,
# weighted median, leave-one-out, Steiger, F-stat — are NOT part of the paper's
# genome-wide MR and are deliberately out of scope here.)
#
# METHOD, step by step
#   1. Candidate instruments = {TRAIT} lead variants (Step-3 output) with
#      GWAS MAF > MAF_MIN and P < GWS_P.
#   2. Greedy LD pruning (verbatim port of the authors' prune_by_LD, from
#      code/MR/cis_MR.R in github.com/ralf-tambets/EstBB-UKBB-metaanalysis):
#        (1) move the lowest-P candidate into the instrument set;
#        (2) drop it and every candidate in LD (r2 > R2) with it;
#        (3) repeat until the candidate pool is empty.
#      r2 is computed ON THE FLY from the 1000G-EUR panel with plink2 — one
#      square r matrix per chromosome (cross-chromosome r2 is 0 by construction),
#      assembled block-diagonally. Leads absent from the panel keep r2 = 0 to all
#      others (survive as independent instruments) — reported as a coverage count.
#   3. Harmonize the pruned instruments against the {OUTCOME} GWAS on rsID:
#      align to the exposure effect allele (flip outcome beta on allele swap),
#      drop strand-ambiguous palindromes (A/T, C/G).
#   4. IVW-MR: mr_input(bx, bxse, by, byse) -> mr_ivw(model = "random").
#
# INPUTS (env vars; self-bridged via `conda run -n variantbench-r` in the rule)
#   TRAIT      exposure/metabolite label            [LDL_C]
#   OUTCOME    outcome label (report/output naming)  [CAD]
#   LEADS      exposure lead-variant table           [results/loci/{TRAIT}.lead.tsv]
#   OUTCOME_GZ outcome GWAS sumstats (gz, harmonized) [data/raw/CAD.GCST90132314.h.tsv.gz]
#   REFDIR     1000G-EUR plink2 panel dir            [data/ref/1kg_phase3_hg38]
#   KEEP       EUR sample keep-file                  [{REFDIR}/1kg_eur.keep]
#   PLINK2     plink2 launcher                        [conda run -n plink2 plink2]
#   OUTDIR     output dir                            [results/mr]
#   GWS_P      genome-wide significance for leads     [5e-8]
#   MAF_MIN    instrument MAF floor                   [0.01]
#   R2         greedy-pruning LD r2 threshold         [0.01]
#   SEED       RNG seed (reproducibility)             [1123]
#
# OUTPUTS
#   {OUTDIR}/{TRAIT}.{OUTCOME}.gwmr.tsv         IVW result (one row)
#   {OUTDIR}/{TRAIT}.{OUTCOME}.instruments.tsv  harmonized per-SNP bx/by table
#   {OUTDIR}/{TRAIT}.{OUTCOME}.scatter.png      per-SNP by vs bx with IVW slope
# =============================================================================
suppressMessages({ library(data.table); library(MendelianRandomization); library(ggplot2) })

env <- function(k, d="") Sys.getenv(k, unset=d)
TRAIT      <- env("TRAIT",   "LDL_C")
OUTCOME    <- env("OUTCOME", "CAD")
LEADS      <- env("LEADS",   sprintf("results/loci/%s.lead.tsv", TRAIT))
OUTCOME_GZ <- env("OUTCOME_GZ", "data/raw/CAD.GCST90132314.h.tsv.gz")
REFDIR     <- env("REFDIR",  "data/ref/1kg_phase3_hg38")
KEEP       <- env("KEEP",    file.path(REFDIR, "1kg_eur.keep"))
PLINK2     <- env("PLINK2",  "conda run -n plink2 plink2")
OUTDIR     <- env("OUTDIR",  "results/mr")
GWS_P      <- as.numeric(env("GWS_P",   "5e-8"))
MAF_MIN    <- as.numeric(env("MAF_MIN", "0.01"))
R2         <- as.numeric(env("R2",      "0.01"))
SEED       <- as.integer(env("SEED",    "1123"))
set.seed(SEED)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("== 06 genome-wide MR | %s -> %s | MAF>%.3g P<%.1e greedy-r2<%.3g ==\n",
            TRAIT, OUTCOME, MAF_MIN, GWS_P, R2))

# ---- 1. candidate instruments: leads, MAF>1%, genome-wide significant --------
leads <- fread(LEADS)
stopifnot(all(c("SNP","chr","pos","ea","oa","maf","beta","se","pval") %in% names(leads)))
cand <- leads[maf > MAF_MIN & pval < GWS_P]
cand[, chr := as.character(chr)]
# Autosomes only: the 1000G-EUR panel is autosomal (chr1-22); the X chromosome
# (coded 23 in the lead table) has a distinct LD/ploidy structure and is
# conventionally excluded from genome-wide MR. Drop it with a logged count so a
# handful of X-linked leads don't silently survive pruning as "independent".
n_pre <- nrow(cand)
cand  <- cand[chr %in% as.character(1:22)]
n_x   <- n_pre - nrow(cand)
# stage tracker: append to a .progress file on disk (visible even when stdout is
# buffered) and flush stdout so the terminal stays live too.
.PROG <- file.path(OUTDIR, sprintf("%s.%s.progress", TRAIT, OUTCOME))
if (file.exists(.PROG)) file.remove(.PROG)
say <- function(msg) {
  cat(msg); flush(stdout())
  cat(sprintf("%s  %s", format(Sys.time(), "%H:%M:%S"), sub("\n$","\n",msg)),
      file = .PROG, append = TRUE)
}
say(sprintf("[1] leads=%d  ->  candidates (MAF>%.3g, P<%.1e)=%d  (dropped %d non-autosomal)\n",
            nrow(leads), MAF_MIN, GWS_P, nrow(cand), n_x))
if (nrow(cand) < 2) stop("fewer than 2 candidate instruments — cannot run IVW")

# ---- 2a. r2 among candidates, per chromosome, from the 1000G-EUR panel --------
# One plink2 --r-unphased square call per chromosome touched by the candidates.
# We point plink2 at the rsID-keyed pvar (chr{N}_hg38_rs.pvar.zst) so --extract
# matches on the leads' rsIDs; r is squared to r2. Cross-chromosome pairs are 0.
ids   <- cand$SNP
n     <- length(ids)
r2mat <- matrix(0, n, n, dimnames = list(ids, ids))
diag(r2mat) <- 1
matched <- character(0)

for (c in sort(unique(cand$chr))) {
  chr_ids <- cand[chr == c, SNP]
  if (length(chr_ids) < 2) { matched <- c(matched, intersect(chr_ids, ids)); next }
  tag <- file.path(OUTDIR, sprintf(".ld_%s_chr%s", TRAIT, c))
  writeLines(chr_ids, paste0(tag, ".ids"))
  pgen <- file.path(REFDIR, sprintf("chr%s_hg38.pgen",       c))
  pvar <- file.path(REFDIR, sprintf("chr%s_hg38_rs.pvar.zst", c))
  psam <- file.path(REFDIR, sprintf("chr%s_hg38.psam",       c))
  if (!file.exists(pgen)) { cat(sprintf("   [warn] panel missing for chr%s — its leads stay independent\n", c)); next }
  cmd <- sprintf(paste("%s --pgen %s --pvar %s --psam %s --keep %s",
                       "--extract %s --r-unphased square ref-based --out %s"),
                 PLINK2, shQuote(pgen), shQuote(pvar), shQuote(psam), shQuote(KEEP),
                 shQuote(paste0(tag, ".ids")), shQuote(tag))
  st <- suppressWarnings(system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE))
  vcor <- paste0(tag, ".unphased.vcor1"); vars <- paste0(tag, ".unphased.vcor1.vars")
  if (st != 0 || !file.exists(vcor) || !file.exists(vars)) {
    cat(sprintf("   [warn] plink2 r failed for chr%s — its leads stay independent\n", c)); next
  }
  v  <- readLines(vars)
  rr <- as.matrix(fread(vcor, header = FALSE))
  dimnames(rr) <- list(v, v)
  rr[is.na(rr)] <- 0
  common <- intersect(v, chr_ids)
  r2mat[common, common] <- rr[common, common]^2          # signed r -> r2
  matched <- c(matched, common)
  file.remove(Sys.glob(paste0(tag, "*")))                 # tidy scratch
}
n_missing <- length(setdiff(ids, matched))
say(sprintf("[2] LD panel coverage: %d/%d candidates matched; %d absent (kept independent)\n",
            length(unique(matched)), n, n_missing))

# plink2 --r emits NaN for a variant that is monomorphic in the 525-EUR panel
# (zero variance -> correlation undefined). Such a NaN lands on that variant's
# DIAGONAL via rr^2, overwriting the 1 -> in prune_by_LD, `which(r2mat[lead,] >
# r2)` then never contains `lead` (NaN > x is NA), so `lead` is never dropped and
# the greedy while-loop spins forever. Force NaN/NA -> 0 and restore diag = 1.
n_bad <- sum(!is.finite(r2mat))
r2mat[!is.finite(r2mat)] <- 0
diag(r2mat) <- 1
if (n_bad > 0) say(sprintf("[2*] cleaned %d non-finite r2 cells (monomorphic-in-panel); diag restored\n", n_bad))

# ---- 2b. greedy LD pruning (verbatim algorithm from authors' prune_by_LD) ----
# sort ascending P, take top, drop it + all r2>R2, repeat until pool empty.
prune_by_LD <- function(cand, r2mat, r2_limit) {
  s <- cand[order(cand$pval)]                              # lowest P first
  keep <- character(0)
  while (nrow(s) > 0) {
    lead <- s$SNP[1]
    keep <- c(keep, lead)
    in_ld <- names(which(r2mat[lead, ] > r2_limit))        # includes self (r2=1)
    s <- s[!(SNP %in% in_ld)]
  }
  cand[SNP %in% keep]
}
say("[2a->] entering greedy prune\n")
instr <- prune_by_LD(cand, r2mat, R2)
say(sprintf("[2b] greedy pruning: %d candidates -> %d independent instruments\n",
            nrow(cand), nrow(instr)))

# ---- 3. harmonize instruments against the outcome GWAS ------------------------
# Match on rsID; align outcome to the EXPOSURE effect allele; drop palindromes.
# The outcome GWAS (~19.5M rows, 1.2 GB BGZF) is tabix-indexed on chr/pos, so we
# fetch ONLY the pruned instruments' loci (a few-KB indexed read) instead of a
# single-threaded full-file gz decompression. Same build (GRCh38) both sides, so
# the position slice returns the identical rows a full read would merge; rsID is
# retained (col 26) so the rsID merge / harmonization below is byte-for-byte the
# same as reading the whole file. (If the .tbi is absent the pipeline builds it:
# `tabix -s 1 -b 2 -e 2 -S 1 <outcome>.h.tsv.gz`.)
read_outcome_at <- function(pos_dt) {
  regs <- sprintf("%s:%d-%d", sub("^chr","",as.character(pos_dt$chr)), pos_dt$pos, pos_dt$pos)
  txt  <- tryCatch(system2("tabix", c(shQuote(OUTCOME_GZ), regs), stdout = TRUE, stderr = FALSE),
                   error = function(e) character(0))
  if (length(txt) == 0) return(data.table())
  d <- fread(text = paste(txt, collapse = "\n"), header = FALSE, sep = "\t")
  # CAD cols: 1 chr | 2 pos | 3 effect_allele | 4 other_allele | 5 beta | 6 se |
  #           7 eaf | 8 p_value | ... | 26 rsid  (BGZF is CRLF -> strip trailing \r)
  d <- d[, .(SNP = sub("\r$","", as.character(V26)),
             o_ea = V3, o_oa = V4, o_beta = as.numeric(V5),
             o_se = as.numeric(V6), o_pval = as.numeric(V8))]
  unique(d[SNP != "" & SNP != "NA"], by = "SNP")
}
say(sprintf("[3->] entering tabix harmonize on %d pruned instruments\n", nrow(instr)))
oc <- read_outcome_at(instr)
h  <- merge(instr[, .(SNP, ea, oa, e_beta = beta, e_se = se)], oc, by = "SNP")
up <- function(x) toupper(as.character(x))
h[, `:=`(ea = up(ea), oa = up(oa), o_ea = up(o_ea), o_oa = up(o_oa))]

# strand-ambiguous palindromes (A/T, C/G) — unresolvable without frequency, drop
comp <- c(A="T", T="A", C="G", G="C")
h <- h[!(o_ea == comp[ea] & o_oa == comp[oa] | ea == oa)]

# allele alignment to exposure effect allele
h <- h[(ea == o_ea & oa == o_oa) | (ea == o_oa & oa == o_ea)]   # matching allele pair only
h[, flip := as.integer(ea == o_oa & oa == o_ea)]                # outcome swapped -> flip sign
h[, by_aligned := ifelse(flip == 1L, -o_beta, o_beta)]
n_dropped <- nrow(instr) - nrow(h)
say(sprintf("[3] harmonized instruments: %d (dropped %d: palindromic / allele-mismatch / absent in outcome)\n",
            nrow(h), n_dropped))
if (nrow(h) < 2) stop("fewer than 2 harmonized instruments — cannot run IVW")

# ---- 4. multiplicative random-effects IVW-MR ---------------------------------
say(sprintf("[4->] entering IVW on %d harmonized instruments\n", nrow(h)))
mri <- mr_input(exposure = TRAIT, outcome = OUTCOME, snps = h$SNP,
                bx = h$e_beta, bxse = h$e_se, by = h$by_aligned, byse = h$o_se)
ivw <- mr_ivw(mri, model = "random")                       # multiplicative random-effects
res <- data.table(
  exposure = TRAIT, outcome = OUTCOME, method = "IVW (mult. random-effects)",
  nsnp     = ivw@SNPs,
  b        = ivw@Estimate,     se = ivw@StdError,
  ci_lo    = ivw@CILower,      ci_up = ivw@CIUpper,
  pval     = ivw@Pvalue,
  or       = exp(ivw@Estimate),
  or_lo    = exp(ivw@CILower), or_up = exp(ivw@CIUpper))
fwrite(res, file.path(OUTDIR, sprintf("%s.%s.gwmr.tsv", TRAIT, OUTCOME)), sep = "\t")
fwrite(h[, .(SNP, ea, oa, bx = e_beta, bxse = e_se, by = by_aligned, byse = o_se,
             flip, o_pval)],
       file.path(OUTDIR, sprintf("%s.%s.instruments.tsv", TRAIT, OUTCOME)), sep = "\t")
say(sprintf("[4] IVW: b=%.4f se=%.4f  OR=%.3f (%.3f, %.3f)  P=%.2e  nsnp=%d\n",
            res$b, res$se, res$or, res$or_lo, res$or_up, res$pval, res$nsnp))

# ---- scatter: per-SNP outcome effect vs exposure effect, IVW slope -----------
pd <- copy(h)
sgn <- sign(pd$e_beta); pd[, `:=`(bx = e_beta * sgn, by = by_aligned * sgn)]  # to right half-plane
say("[5->] entering scatter plot\n")
p <- ggplot(pd, aes(bx, by)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_errorbar(aes(ymin = by - o_se, ymax = by + o_se), colour = "grey60", width = 0, linewidth = 0.3) +
  geom_errorbarh(aes(xmin = bx - e_se, xmax = bx + e_se), colour = "grey60", height = 0, linewidth = 0.3) +
  geom_abline(slope = res$b, intercept = 0, colour = "#c0392b", linewidth = 0.8) +
  geom_point(size = 1.6, colour = "#2c3e50") +
  labs(x = sprintf("Per-variant effect on %s (SD)", TRAIT),
       y = sprintf("Per-variant effect on %s (log-OR)", OUTCOME),
       title = sprintf("Genome-wide IVW-MR: %s -> %s", TRAIT, OUTCOME),
       subtitle = sprintf("%d instruments  |  OR=%.2f (%.2f-%.2f)  P=%.1e",
                          res$nsnp, res$or, res$or_lo, res$or_up, res$pval)) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 10))
ggsave(file.path(OUTDIR, sprintf("%s.%s.scatter.png", TRAIT, OUTCOME)),
       plot = p, width = 5.2, height = 5.0, dpi = 150, bg = "white")
say("[done] 06 genome-wide MR\n")
