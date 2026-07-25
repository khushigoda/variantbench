#!/usr/bin/env Rscript
# Step 0 — INSPECTION (standalone; NOT a Snakemake rule, produces no analysis data).
# It reads only the HEADER + a small sample of rows from each raw EstBB / UKBB_EUR file
# (kilobytes, never a full-file load) to answer the three structural questions that
# determine how 01 and 02 are written. Output is a short findings note we commit.
# meta_EUR is NOT inspected — it is Step 1's validation target and is checked against
# the paper's reported numbers, not combined, so its schema doesn't shape 01/02 code.
#
#   Run from the repo root:  Rscript workflow/scripts/00_inspect.R
#
# Why this is out of the DAG: it produces nothing that 01/02 consume as input. Its
# findings travel as decisions baked into 01/02's code, not as a data file. See
# ANALYSIS.md / the Snakefile comment block.
#
# The three questions:
#   Q1  Is the GWAS-SSF column schema identical across all 8 files (4 traits x 2 cohorts)?
#       -> if yes, 01/02 use one uniform column-selective read.
#   Q2  What does `variant_id` hold — rsID (rs...) or positional (chr_pos_ref_alt)?
#       -> LDSC munge merges to w_hm3.snplist BY rsID. Positional => 02 must attach
#          rsIDs before munge; rsID => munge merges directly.
#   Q3  Is chromosome coding consistent (numeric vs 'chr' prefix; X vs 23)?
#       -> so 01's join keys (chr:pos:ea:oa) align across cohorts.
suppressMessages({library(data.table); library(yaml)})

cfg     <- yaml::read_yaml("config/config.yaml")
metabs  <- names(cfg$metabolites)
cohorts <- cfg$cohorts                            # EstBB, UKBB_EUR only — the two cohorts
                                                  # Step 1 combines. meta_EUR is a validation
                                                  # target (checked vs the paper), not inspected.
NSAMPLE <- 2000L                                  # rows sampled per file (never full read)
dir.create("results/inspect", recursive = TRUE, showWarnings = FALSE)

raw_path <- function(m, co) {
  acc <- cfg$metabolites[[m]]$accession[[co]]
  file.path("data/raw", sprintf("%s_%s.%s.tsv.gz", m, co, acc))
}

rows <- list()
for (m in metabs) for (co in cohorts) {
  f <- raw_path(m, co)
  if (!file.exists(f)) {
    rows[[length(rows) + 1]] <- data.table(
      metab = m, cohort = co, file = basename(f), exists = FALSE,
      ncol = NA_integer_, header = NA_character_,
      id_example = NA_character_, id_type = NA_character_, chr_example = NA_character_)
    next
  }
  s <- tryCatch(fread(f, nrows = NSAMPLE, showProgress = FALSE), error = function(e) NULL)
  if (is.null(s)) next
  idcol  <- intersect(c("variant_id", "rsid", "rs_id", "ID"), names(s))
  idex   <- if (length(idcol)) as.character(s[[idcol[1]]][1]) else NA_character_
  idtype <- if (is.na(idex)) NA_character_
            else if (grepl("^rs[0-9]+", idex)) "rsID" else "positional/other"
  chrcol <- intersect(c("chromosome", "chr", "CHR", "CHROM"), names(s))
  chrex  <- if (length(chrcol)) as.character(s[[chrcol[1]]][1]) else NA_character_
  rows[[length(rows) + 1]] <- data.table(
    metab = m, cohort = co, file = basename(f), exists = TRUE,
    ncol = ncol(s), header = paste(names(s), collapse = ","),
    id_example = idex, id_type = idtype, chr_example = chrex)
}
rep <- rbindlist(rows, fill = TRUE)
fwrite(rep, "results/inspect/schema_report.tsv", sep = "\t")

# ---- aggregate the three findings ----
present <- rep[exists == TRUE]
schemas <- unique(present$header)
q1 <- if (nrow(present) == 0) "no files present yet — re-run once downloads finish"
      else if (length(schemas) == 1) "IDENTICAL across all present files"
      else sprintf("DIFFERS: %d distinct header layouts (see schema_report.tsv)", length(schemas))
q2 <- paste(unique(na.omit(present$id_type)), collapse = ", ")
q3 <- paste(unique(na.omit(present$chr_example)), collapse = ", ")
q2_action <- if (grepl("positional", q2)) "02 MUST attach rsIDs before LDSC munge"
             else if (identical(q2, "rsID")) "SNP col is rsID; munge merges directly"
             else "unknown — re-check once files land"

md <- c(
  "# Step 0 — schema inspection findings",
  sprintf("_Generated %s. Sampled %d rows/file. No full-file reads; no analysis data written._",
          Sys.Date(), NSAMPLE),
  "",
  sprintf("- **Files present:** %d / %d expected (4 traits x {EstBB, UKBB_EUR}).",
          nrow(present), nrow(rep)),
  sprintf("- **Q1 — schema uniform?** %s", q1),
  sprintf("- **Q2 — variant_id type:** %s  ->  %s", q2, q2_action),
  sprintf("- **Q3 — chr coding seen:** %s", q3),
  "",
  "Per-file detail: `results/inspect/schema_report.tsv`.")
writeLines(md, "results/inspect/00_findings.md")
cat(paste(md, collapse = "\n"), "\n")
