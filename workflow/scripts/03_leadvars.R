#!/usr/bin/env Rscript
# =====================================================================================
# Step 3 - Identification of lead variants (genome-wide-significant loci) + Manhattan.
#
#   snakemake@input$meta  : results/meta/<metab>.meta.tsv.gz  (Step-1 output, ~104M rows)
#   snakemake@output$lead : results/loci/<metab>.lead.tsv
#   snakemake@output$plot : results/loci/<metab>.manhattan.png
#
# WHY THIS IS STREAMED, NOT fread()'d:
#   The meta table is ~104M variants / ~5 GB gzipped per trait. Loading it whole into R
#   on an 8 GB machine OOMs. Instead a single `zcat | awk` pass pre-filters the file:
#   only the (few hundred) significant hits and a thinned point-set for the plot ever
#   cross into R. Memory stays flat regardless of trait/genome size.
#
# SIGNIFICANCE - DUAL THRESHOLD BY ALLELE FREQUENCY (paper's method):
#   Common / low-freq variants (MAF > maf_rare, default 0.1%) use the standard
#     genome-wide line P < 5e-8.
#   Rare variants (MAF < maf_rare) use a stricter Bonferroni line P < 0.05/8e7 = 6.25e-10,
#     because ~80M rare tests make 5e-8 too lenient (≈4 expected false positives).
#   We threshold on neg_log10p (computed in log space upstream) rather than pval, so
#   ultra-strong loci (pval underflows to 0) are handled correctly.
#   eaf-unknown variants (both cohorts missing frequency) fall back to the standard line
#     and are flagged maf_class = "unknown" so the report can caveat them.
#
# LOCUS DEFINITION - DISTANCE-BASED CLUMPING (no LD panel):
#   The brief's Step 3 is summary-statistics-only and provides no genotype panel, so we
#   cannot compute pairwise r2. We use the standard distance rule the paper's window is
#   built on: greedily take the most-significant unclaimed variant as a lead, claim a
#   locus_window-wide locus around it (default 2 Mb, ±1 Mb), remove everything inside,
#   iterate. The paper's extra r2>=0.05 locus-merge and r2>=0.8 cross-trait clustering
#   need a EUR genotype panel and are deferred to the LD-aware follow-on (Step 4 sourcing).
# =====================================================================================
suppressMessages({library(data.table); library(ggplot2)})
set.seed(snakemake@params$seed)                 # reproducibility (thinning is deterministic anyway)

meta_gz  <- snakemake@input$meta
gws      <- as.numeric(snakemake@params$gws_p)        # 5e-8
rare_gws <- as.numeric(snakemake@params$rare_gws_p)   # 6.25e-10
maf_rare <- as.numeric(snakemake@params$maf_rare)     # 0.001
win      <- as.numeric(snakemake@params$locus_window) # 2e6
half     <- win / 2
trait    <- snakemake@wildcards$metabolite

thr_c   <- -log10(gws)        # 7.301  common / unknown
thr_r   <- -log10(rare_gws)   # 9.204  rare
MANH_FLOOR <- 3               # keep every point with -log10p >= 3 (p<1e-3): all real signal
THIN       <- 500             # plus 1-in-THIN of the rest as background "grass"

# ---- meta column order (from 01_meta.R, authoritative) ----
# 1 SNP 2 chr 3 pos 4 ea 5 oa 6 eaf 7 beta 8 se 9 z 10 pval 11 neg_log10p 12 n
# 13 direction 14 het_q 15 het_p 16 n_cohorts

hits_f <- tempfile(fileext = ".tsv")             # significant rows, full 16 cols
pts_f  <- tempfile(fileext = ".tsv")             # thinned chr/pos/nlp/SNP for the plot

# ---- ONE streaming pass: dual-threshold filter + Manhattan thinning ----
awk_prog <- '
BEGIN{ FS=OFS="\t" }
NR==1 { print $0 > hits; print "chr","pos","nlp","SNP" > pts; next }
{
  eaf=$6; nlp=$11+0
  if (eaf=="NA" || eaf=="") { maf=-1 }
  else { e=eaf+0; maf=(e<0.5? e : 1-e) }
  if      (maf<0)          thr=thr_c      # unknown -> standard line (flagged in R)
  else if (maf>=maf_rare)  thr=thr_c      # common
  else                     thr=thr_r      # rare
  if (nlp > thr) print $0 > hits
  if (nlp >= manh_floor)   { print $2,$3,$11,$1 > pts }
  else { c++; if (c % thin == 0) print $2,$3,$11,$1 > pts }
}'

cmd <- sprintf(
  "zcat < %s | awk -v hits=%s -v pts=%s -v thr_c=%f -v thr_r=%f -v maf_rare=%f -v manh_floor=%d -v thin=%d %s",
  shQuote(meta_gz), shQuote(hits_f), shQuote(pts_f),
  thr_c, thr_r, maf_rare, MANH_FLOOR, THIN, shQuote(awk_prog))
# macOS `zcat` refuses .gz without a flag; `zcat <` (stdin) is portable. gzip -dc fallback:
cmd <- sub("^zcat < ", "gzip -dc ", cmd)
message(sprintf("[%s] streaming filter over %s ...", trait, basename(meta_gz)))
st <- system(cmd)
if (st != 0) stop(sprintf("awk stream failed (exit %d) for %s", st, meta_gz))

hits <- fread(hits_f)
pts  <- fread(pts_f)
unlink(c(hits_f, pts_f))
message(sprintf("[%s] %d significant variants, %d thinned plot points",
                trait, nrow(hits), nrow(pts)))

# ---- classify each significant hit (for reporting / MCQ) ----
if (nrow(hits)) {
  hits[, maf := fifelse(is.na(eaf), NA_real_, pmin(eaf, 1 - eaf))]
  hits[, maf_class := fifelse(is.na(maf), "unknown",
                       fifelse(maf >= maf_rare, "common", "rare"))]
}

# ---- distance-based greedy clumping -> lead variants ----
clump_chr <- function(d, half) {
  setorder(d, -neg_log10p)                       # most significant first
  leadpos <- numeric(0); keep <- logical(nrow(d))
  for (i in seq_len(nrow(d))) {
    p <- d$pos[i]
    if (length(leadpos) == 0 || all(abs(leadpos - p) > half)) {
      keep[i] <- TRUE; leadpos <- c(leadpos, p)
    }
  }
  d[keep]
}
if (nrow(hits)) {
  leads <- rbindlist(lapply(split(hits, by = "chr"), clump_chr, half = half))
  setorder(leads, chr, pos)
  # locus bounds + how many significant variants fall in each lead's window
  leads[, `:=`(locus_start = pmax(0, pos - half), locus_end = pos + half)]
  leads[, n_sig_in_locus := {
    cc <- chr; ps <- pos
    vapply(seq_len(.N), function(i)
      hits[chr == cc[i] & abs(pos - ps[i]) <= half, .N], integer(1))
  }]
  lead_out <- leads[, .(SNP, chr, pos, ea, oa, eaf, maf, maf_class,
                        beta, se, z, pval, neg_log10p, n, direction, het_p,
                        n_cohorts, locus_start, locus_end, n_sig_in_locus)]
} else {
  lead_out <- data.table(SNP=character(), chr=integer(), pos=numeric(), ea=character(),
                         oa=character(), eaf=numeric(), maf=numeric(), maf_class=character(),
                         beta=numeric(), se=numeric(), z=numeric(), pval=numeric(),
                         neg_log10p=numeric(), n=numeric(), direction=character(),
                         het_p=numeric(), n_cohorts=integer(), locus_start=numeric(),
                         locus_end=numeric(), n_sig_in_locus=integer())
}
fwrite(lead_out, snakemake@output$lead, sep = "\t")
message(sprintf("[%s] %d independent lead loci (common=%d rare=%d unknown=%d)",
                trait, nrow(lead_out),
                sum(lead_out$maf_class=="common"), sum(lead_out$maf_class=="rare"),
                sum(lead_out$maf_class=="unknown")))

# ---- Manhattan plot: cumulative genome position, alternating chr bands, leads marked ----
pts <- pts[is.finite(nlp)]
setorder(pts, chr, pos)
chr_mx  <- pts[, .(mx = max(pos)), by = chr][order(chr)]
chr_mx[, offset := cumsum(as.numeric(shift(mx, fill = 0)))]   # as.numeric: genome cumpos ~3.1e9 overflows 32-bit int
pts <- merge(pts, chr_mx[, .(chr, offset)], by = "chr")
pts[, xpos := pos + offset]
pts[, band := factor(chr %% 2)]
axis_df <- merge(chr_mx, pts[, .(center = median(xpos)), by = chr], by = "chr")

pl <- ggplot(pts, aes(xpos, nlp, colour = band)) +
  geom_point(size = 0.35, alpha = 0.7) +
  scale_colour_manual(values = c("0" = "grey65", "1" = "grey35"), guide = "none") +
  geom_hline(yintercept = thr_c, linetype = 2, colour = "red") +
  geom_hline(yintercept = thr_r, linetype = 3, colour = "purple") +
  scale_x_continuous(breaks = axis_df$center,
                     labels = ifelse(axis_df$chr == 23, "X", as.character(axis_df$chr))) +
  labs(x = "Chromosome", y = expression(-log[10](italic(P))),
       title = sprintf("%s - meta-analysis (EstBB + UKBB_EUR)", trait),
       subtitle = sprintf("%d lead loci; red = 5e-8, purple = rare 6.25e-10",
                          nrow(lead_out))) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

if (nrow(lead_out)) {
  ld <- merge(lead_out[, .(chr, pos, SNP, nlp = neg_log10p)],
              chr_mx[, .(chr, offset)], by = "chr")
  ld[, xpos := pos + offset]
  pl <- pl +
    geom_point(data = ld, aes(xpos, nlp), colour = "firebrick", size = 1.6, inherit.aes = FALSE)
  topn <- head(ld[order(-nlp)], 8)               # label the 8 strongest loci
  pl <- pl + geom_text(data = topn, aes(xpos, nlp, label = SNP),
                       inherit.aes = FALSE, size = 2.4, vjust = -0.6, check_overlap = TRUE)
}
ggsave(snakemake@output$plot, pl, width = 11, height = 4.5, dpi = 130)
message(sprintf("[%s] wrote %s + %s", trait,
                basename(snakemake@output$lead), basename(snakemake@output$plot)))

# ---- MAF vs |effect| scatter (analogue of Karjalainen et al. Supp. Fig. 3A) --------
#   Shows the inverse allele-frequency/effect-size relationship: rare lead variants
#   carry large effects, common leads small ones. Point size = -log10 P; colour = the
#   frequency class used for the dual significance threshold. Only leads with a KNOWN
#   MAF class enter (unknown-frequency leads are excluded, not silently binned).
if (nrow(lead_out)) {
  me <- lead_out[maf_class %in% c("common","rare") & is.finite(maf) & is.finite(beta) & maf > 0]
  if (nrow(me)) {
    me[, absbeta := abs(beta)]
    me[, freq := factor(ifelse(maf_class=="rare","rare (MAF<0.1%)","common (MAF>0.1%)"),
                        levels=c("common (MAF>0.1%)","rare (MAF<0.1%)"))]
    # label PCSK9's R46L drug-target lead if it is itself a sentinel (it is, for LDL_C)
    lab_me <- me[SNP == "rs11591147"]
    if (nrow(lab_me)) lab_me[, gene := "PCSK9"]

    p_me <- ggplot(me, aes(maf, absbeta)) +
      geom_point(aes(size = neg_log10p, colour = freq), alpha = 0.6) +
      scale_x_log10(breaks = c(1e-4,1e-3,1e-2,1e-1,0.5),
                    labels = c("0.01%","0.1%","1%","10%","50%"),
                    expand = expansion(mult = c(0.08, 0.05))) +
      scale_size_continuous(range = c(0.6, 5), name = expression(-log[10](italic(P)))) +
      scale_colour_manual(values = c("common (MAF>0.1%)"="#4C72B0","rare (MAF<0.1%)"="#DD8452"),
                          name = NULL) +
      labs(x = "Lead-variant minor allele frequency (log scale)",
           y = expression("|effect size|  |"*beta*"|"),
           title = sprintf("Rare %s lead variants carry larger effect sizes", trait),
           subtitle = sprintf("meta_EUR (EstBB + UKBB_EUR); %d independent leads", nrow(me))) +
      theme_bw(base_size = 11) +
      theme(panel.grid.minor = element_blank(), legend.position = "right",
            plot.title = element_text(face = "plain"))
    if (nrow(lab_me)) {
      p_me <- p_me +
        geom_point(data = lab_me, aes(maf, absbeta), shape = 21, size = 2.6,
                   colour = "black", fill = NA, stroke = 0.5, inherit.aes = FALSE) +
        geom_text(data = lab_me, aes(maf, absbeta, label = gene), size = 3,
                  fontface = "italic", vjust = -0.9, inherit.aes = FALSE) +
        labs(subtitle = sprintf("meta_EUR (EstBB + UKBB_EUR); %d independent leads; PCSK9 = fine-mapped coding lead", nrow(me)))
    }
    ggsave(snakemake@output$maf_effect, p_me, width = 8, height = 5, dpi = 140)
    message(sprintf("[%s] wrote %s", trait, basename(snakemake@output$maf_effect)))
  }
}
