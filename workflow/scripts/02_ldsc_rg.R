#!/usr/bin/env Rscript
# Step 2 (rg + parse + plots) — LDSC wrapper, paper replication over 3 datasets x 4 traits.
# munge (-> .sumstats.gz) and --h2 run as shell rules; this R step runs the two rg analyses
# from the paper, parses the h2 + rg LOG files into tidy tables, and renders the figures.
# All LDSC calls go through the `ldsc` conda env.
#
#   Analysis A — "between biobanks for each trait": for each trait, rg the 3 datasets
#                (EstBB, UKBB_EUR, meta_EUR) against each other -> C(3,2)=3 x 4 = 12 rg.
#   Analysis B — "between all traits in three datasets": within each dataset, rg the 4
#                traits against each other -> C(4,2)=6 x 3 = 18 rg.
#
# Munged files are named results/ldsc/<trait>.<dataset>.sumstats.gz, so (trait, dataset)
# is recovered from the filename. LDSC's `--rg f1,f2,...,fk` correlates f1 vs each of the
# rest (a "staircase"), so to cover every C(k,2) pair once we lead each element i over the
# tail i..k. One short LDSC call per lead; logs are parsed and combined.
suppressMessages({library(data.table); library(ggplot2)})
set.seed(snakemake@params$seed)

sumstats <- unlist(snakemake@input$sumstats)
h2logs   <- unlist(snakemake@input$h2)
ld       <- snakemake@params$ld
ldsc_env <- snakemake@params$ldsc_env
metabs   <- snakemake@params$metabs      # trait order  (config: metabolites)
datasets <- snakemake@params$datasets    # dataset order (config: ldsc_datasets)
dir.create("results/ldsc", showWarnings = FALSE, recursive = TRUE)

# key = "<trait>.<dataset>"; map key -> munged file path
keys  <- sub("\\.sumstats\\.gz$", "", basename(sumstats))
paths <- setNames(sumstats, keys)
key_of <- function(trait, ds) paste(trait, ds, sep = ".")

# ---- helper: run one staircase of rg over an ordered vector of keys ----
# Returns the path of the produced .log (one call: lead = first key vs the rest).
run_rg <- function(ordered_keys, out_prefix) {
  files <- paths[ordered_keys]
  system2("conda", c("run", "-n", ldsc_env, "python", path.expand("~/ldsc/ldsc.py"),
                     "--rg", paste(files, collapse = ","),
                     "--ref-ld-chr", ld, "--w-ld-chr", ld,
                     "--out", out_prefix))
  paste0(out_prefix, ".log")
}

# ---- rg summary-block parser (whitespace-aligned; read.table splits on any run) ----
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
  tryCatch(as.data.table(read.table(text = paste(blk, collapse = "\n"),
                                     header = TRUE, stringsAsFactors = FALSE)),
           error = function(e) data.table())
}
# turn a parsed rg block's p1/p2 file paths back into "<trait>.<dataset>" keys
clean_keys <- function(dt) {
  if (!nrow(dt)) return(dt)
  dt[, p1 := sub("\\.sumstats\\.gz$", "", basename(p1))]
  dt[, p2 := sub("\\.sumstats\\.gz$", "", basename(p2))]
  dt[]
}

# ---- Analysis A: between biobanks, per trait ----
# For each trait, staircase over datasets (i..k) so every dataset pair is covered once.
rgA <- list()
for (t in metabs) {
  ks <- key_of(t, datasets)                       # e.g. BCAA.EstBB, BCAA.UKBB_EUR, BCAA.meta_EUR
  ks <- ks[ks %in% keys]
  for (i in seq_len(length(ks) - 1L)) {
    lg <- run_rg(ks[i:length(ks)],
                 sprintf("results/ldsc/rgA_%s_%s", t, sub(".*\\.", "", ks[i])))
    d <- clean_keys(parse_rg(lg))
    if (nrow(d)) { d[, trait := t]; rgA[[length(rgA) + 1L]] <- d }
  }
}
rgA <- rbindlist(rgA, fill = TRUE)
if (nrow(rgA)) {
  rgA[, ds1 := sub(".*\\.", "", p1)][, ds2 := sub(".*\\.", "", p2)]
  rgA <- unique(rgA, by = c("trait", "ds1", "ds2"))
}
fwrite(rgA, snakemake@output$rg_bb, sep = "\t")

# ---- Analysis B: between traits, per dataset ----
rgB <- list()
for (ds in datasets) {
  ks <- key_of(metabs, ds)
  ks <- ks[ks %in% keys]
  for (i in seq_len(length(ks) - 1L)) {
    lg <- run_rg(ks[i:length(ks)],
                 sprintf("results/ldsc/rgB_%s_%s", ds, sub("\\..*", "", ks[i])))
    d <- clean_keys(parse_rg(lg))
    if (nrow(d)) { d[, dataset := ds]; rgB[[length(rgB) + 1L]] <- d }
  }
}
rgB <- rbindlist(rgB, fill = TRUE)
if (nrow(rgB)) {
  rgB[, tr1 := sub("\\..*", "", p1)][, tr2 := sub("\\..*", "", p2)]
  rgB <- unique(rgB, by = c("dataset", "tr1", "tr2"))
}
fwrite(rgB, snakemake@output$rg_tr, sep = "\t")

# ---- parse h2 from each <trait>.<dataset>.h2.log ----
parse_h2 <- function(f) {
  key <- sub("\\.h2\\.log$", "", basename(f))
  trait <- sub("\\..*", "", key); ds <- sub(".*?\\.", "", key)
  if (!file.exists(f)) return(data.table(trait = trait, dataset = ds,
                                          h2 = NA_real_, h2_se = NA_real_, intercept = NA_real_))
  L <- readLines(f, warn = FALSE)
  hit <- grep("Total Observed scale h2:", L, value = TRUE)
  int <- grep("Intercept:", L, value = TRUE)
  h2 <- se <- icpt <- NA_real_
  if (length(hit)) {
    mm <- regmatches(hit[1], regexec("h2:\\s*(-?[0-9.eE+-]+)\\s*\\(([0-9.eE+-]+)\\)", hit[1]))[[1]]
    if (length(mm) == 3) { h2 <- as.numeric(mm[2]); se <- as.numeric(mm[3]) }
  }
  # Anchor on "Intercept:" then capture the number. A naive [0-9.eE+-]+ match fails here:
  # the class contains 'e'/'E' (sci-notation) and would grab the 'e' in the word "Intercept".
  if (length(int)) {
    im <- regmatches(int[1], regexec("Intercept:\\s*(-?[0-9.eE+-]+)", int[1]))[[1]]
    if (length(im) == 2) icpt <- as.numeric(im[2])
  }
  data.table(trait = trait, dataset = ds, h2 = h2, h2_se = se, intercept = icpt)
}
h2tab <- rbindlist(lapply(h2logs, parse_h2), fill = TRUE)
h2tab[, `:=`(trait = factor(trait, metabs), dataset = factor(dataset, datasets))]
fwrite(h2tab, snakemake@output$h2, sep = "\t")

# ---- figures ----
# h2: bar +/- SE, trait on x, faceted by dataset
ph2 <- ggplot(h2tab, aes(trait, h2, fill = dataset)) +
  geom_col(position = position_dodge(width = .8), width = .7) +
  geom_errorbar(aes(ymin = h2 - h2_se, ymax = h2 + h2_se),
                position = position_dodge(width = .8), width = .25) +
  facet_wrap(~ dataset) + theme_bw() + guides(fill = "none") +
  labs(x = NULL, y = expression(paste("SNP ", italic(h)^2)),
       title = "SNP heritability by dataset (LDSC)")
ggsave(snakemake@output$h2_plot, ph2, width = 7, height = 4, dpi = 120)

# symmetric-heatmap helper: rows have (a, b, rg); build full grid with diagonal = 1
heat_facets <- function(dt, a, b, rgcol = "rg", levs, facet, title) {
  if (!nrow(dt)) return(ggplot() + theme_void() + labs(title = paste(title, "— no results")))
  d    <- dt[, .(facet = get(facet), a = get(a), b = get(b), rg = get(rgcol))]
  diag <- CJ(facet = unique(d$facet), a = levs)[, .(facet, a, b = a, rg = 1)]
  dd   <- rbind(d, d[, .(facet, a = b, b = a, rg)], diag)   # off-diag both ways + diagonal
  dd[, `:=`(a = factor(a, levs), b = factor(b, levs))]
  ggplot(dd, aes(a, b, fill = rg)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = ifelse(is.na(rg), "", sprintf("%.2f", rg))), size = 2.8) +
    scale_fill_gradient2(limits = c(-1, 1), low = "#b2182b", mid = "white", high = "#2166ac") +
    facet_wrap(~ facet) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = NULL, title = title)
}
pA <- heat_facets(rgA, "ds1", "ds2", levs = datasets, facet = "trait",
                  title = "Analysis A — rg between biobanks (per trait)")
ggsave(snakemake@output$rg_bb_plot, pA, width = 8, height = 6, dpi = 120)
pB <- heat_facets(rgB, "tr1", "tr2", levs = metabs, facet = "dataset",
                  title = "Analysis B — rg between traits (per dataset)")
ggsave(snakemake@output$rg_tr_plot, pB, width = 9, height = 4, dpi = 120)

cat(sprintf("LDSC: h2 %d rows; rg_between_biobanks %d; rg_between_traits %d\n",
            nrow(h2tab), nrow(rgA), nrow(rgB)))
