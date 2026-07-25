#!/usr/bin/env Rscript
# Step 2 (rg + parse + plots) — LDSC wrapper.
# munge (raw -> .sumstats.gz) and --h2 run as shell rules in the Snakefile; this R step
# runs --rg across all trait pairs, then parses the h2 and rg LOG files into tidy tables
# and renders the two report figures. All LDSC calls go through the `ldsc` conda env.
#
#   input$sumstats : results/ldsc/<metab>.sumstats.gz  (all traits)
#   input$h2       : results/ldsc/<metab>.h2.log       (all traits)
#   output$rg      : results/ldsc/rg.tsv     output$h2 : results/ldsc/h2.tsv
#   output$rg_plot : results/ldsc/rg_heatmap.png       output$h2_plot : h2_barplot.png
suppressMessages({library(data.table); library(ggplot2)})
set.seed(snakemake@params$seed)

sumstats <- unlist(snakemake@input$sumstats)
h2logs   <- unlist(snakemake@input$h2)
ld       <- snakemake@params$ld
ldsc_env <- snakemake@params$ldsc_env
metabs   <- sub("\\.sumstats\\.gz$", "", basename(sumstats))
dir.create("results/ldsc", showWarnings = FALSE, recursive = TRUE)

# ---- 1. run --rg for ALL trait pairs (staircase) ----
# LDSC's `--rg f1,f2,...,fk` computes rg of the FIRST file against each of the rest
# (f1-vs-f2, f1-vs-f3, ...), NOT all pairwise combinations. To cover every C(n,2) pair
# exactly once, lead each trait i over the tail (i+1..n): run i gives rg(i, j) for j>i.
# n-1 short LDSC calls; their logs are concatenated and parsed together below.
rg_prefix <- "results/ldsc/rg"
rg_logs   <- character(0)
n <- length(sumstats)
for (i in seq_len(n - 1L)) {
  out_i <- sprintf("%s_%s", rg_prefix, metabs[i])
  system2("conda", c("run", "-n", ldsc_env, "python", path.expand("~/ldsc/ldsc.py"),
                     "--rg", paste(sumstats[i:n], collapse = ","),
                     "--ref-ld-chr", ld, "--w-ld-chr", ld,
                     "--out", out_i))
  rg_logs <- c(rg_logs, paste0(out_i, ".log"))
}

# ---- 2. parse h2 from each *.h2.log: "Total Observed scale h2: X (SE)" ----
parse_h2 <- function(f) {
  if (!file.exists(f)) return(data.table(metab = sub("\\.h2\\.log$", "", basename(f)),
                                          h2 = NA_real_, h2_se = NA_real_, intercept = NA_real_))
  L <- readLines(f, warn = FALSE)
  hit <- grep("Total Observed scale h2:", L, value = TRUE)
  int <- grep("Intercept:", L, value = TRUE)
  num <- function(s) as.numeric(regmatches(s, regexpr("-?[0-9.]+", s)))
  h2 <- se <- icpt <- NA_real_
  if (length(hit)) {
    mm <- regmatches(hit[1], regexec("h2:\\s*(-?[0-9.eE+-]+)\\s*\\(([0-9.eE+-]+)\\)", hit[1]))[[1]]
    if (length(mm) == 3) { h2 <- as.numeric(mm[2]); se <- as.numeric(mm[3]) }
  }
  if (length(int)) icpt <- num(int[1])
  data.table(metab = sub("\\.h2\\.log$", "", basename(f)), h2 = h2, h2_se = se, intercept = icpt)
}
h2tab <- rbindlist(lapply(h2logs, parse_h2), fill = TRUE)
fwrite(h2tab, snakemake@output$h2, sep = "\t")

# ---- 3. parse the rg summary block from each rg log, then combine ----
# The block is whitespace-ALIGNED (variable-width columns), so read.table (splits on any
# run of whitespace) parses it robustly where fread (single-delimiter) would mis-split.
parse_rg <- function(f) {
  if (!file.exists(f)) return(data.table())
  L <- readLines(f, warn = FALSE)
  i <- grep("Summary of Genetic Correlation Results", L)
  if (!length(i)) return(data.table())
  blk <- L[(i + 1):length(L)]
  stop_at <- which(!grepl("[^[:space:]]", blk))     # first fully-blank line ends the table
  if (length(stop_at)) blk <- blk[seq_len(stop_at[1] - 1L)]
  blk <- blk[grepl("[^[:space:]]", blk)]
  if (length(blk) < 2) return(data.table())          # header + >=1 row required
  dt <- tryCatch(as.data.table(read.table(text = paste(blk, collapse = "\n"),
                                           header = TRUE, stringsAsFactors = FALSE)),
                 error = function(e) data.table())
  dt
}
rgtab <- rbindlist(lapply(rg_logs, parse_rg), fill = TRUE)
if (nrow(rgtab)) {
  rgtab[, p1 := sub("\\.sumstats\\.gz$", "", basename(p1))]
  rgtab[, p2 := sub("\\.sumstats\\.gz$", "", basename(p2))]
  rgtab <- unique(rgtab, by = c("p1", "p2"))         # de-dup any overlap across logs
}
fwrite(rgtab, snakemake@output$rg, sep = "\t")

# ---- 4. figures ----
# h2 bar +/- SE
h2tab[, metab := factor(metab, levels = metabs)]
ph2 <- ggplot(h2tab, aes(metab, h2)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = h2 - h2_se, ymax = h2 + h2_se), width = .2) +
  theme_bw() + labs(x = NULL, y = expression(paste("SNP ", italic(h)^2)),
                    title = "SNP heritability (LDSC)")
ggsave(snakemake@output$h2_plot, ph2, width = 5, height = 4, dpi = 120)

# rg heatmap (symmetric; diagonal = 1)
if (nrow(rgtab)) {
  hm <- rbind(rgtab[, .(p1, p2, rg)],
              rgtab[, .(p1 = p2, p2 = p1, rg)],
              data.table(p1 = metabs, p2 = metabs, rg = 1))
  hm[, `:=`(p1 = factor(p1, metabs), p2 = factor(p2, metabs))]
  prg <- ggplot(hm, aes(p1, p2, fill = rg)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.2f", rg)), size = 3) +
    scale_fill_gradient2(limits = c(-1, 1), low = "#b2182b", mid = "white", high = "#2166ac") +
    theme_bw() + labs(x = NULL, y = NULL, title = "Genetic correlation (LDSC rg)")
} else {
  prg <- ggplot() + theme_void() + labs(title = "rg — no results parsed")
}
ggsave(snakemake@output$rg_plot, prg, width = 5, height = 4.5, dpi = 120)
cat(sprintf("LDSC: parsed h2 for %d traits, rg for %d pairs\n", nrow(h2tab), nrow(rgtab)))
