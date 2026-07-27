#!/usr/bin/env Rscript
# =====================================================================
# Step 8 render harness — self-bridge for `rule report`.
#
# WHY THIS EXISTS: the report is the only pipeline step whose toolchain
# (rmarkdown + DT + pandoc) lives solely in the `variantbench-r` conda
# env, while `snakemake` itself runs in `variantbench`. The old rule used
# `conda:` + `script:`, which only bridges under `--use-conda` — and that
# reintroduces a hanging conda solve on this machine. So `rule report`
# now self-bridges exactly like the coloc/MR rules: its `run:` block
# serialises the rule's REAL resolved inputs to JSON, then calls this
# script via `conda run -n variantbench-r Rscript ...`.
#
# This reconstructs the same S4 `snakemake` object the `script:` directive
# would have injected (input list + params$seed), so 08_report.Rmd is used
# verbatim — no drift between the DAG's paths and the render's paths,
# because both come from the one `input:` block in the Snakefile.
#
# Usage:  Rscript 08_render_report.R <inputs.json> <output.html>
# =====================================================================
suppressMessages({library(rmarkdown); library(methods); library(jsonlite)})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
json_path <- args[[1]]
out_html  <- args[[2]]

payload <- fromJSON(json_path, simplifyVector = FALSE)

# Coerce each input entry to a plain character vector (singletons -> length-1,
# expand() lists -> length-n) — matches how snakemake@input$key behaves.
inp <- lapply(payload$input,  function(x) as.character(unlist(x, use.names = FALSE)))
prm <- lapply(payload$params, function(x) unlist(x, use.names = FALSE))

Snakemake <- setClass("Snakemake",
  representation(input = "list", output = "list",
                 params = "list", wildcards = "list"))
snakemake <- Snakemake(input = inp, output = list(out_html),
                       params = prm, wildcards = list())
assign("snakemake", snakemake, envir = .GlobalEnv)

rmarkdown::render("workflow/scripts/08_report.Rmd",
  output_file   = basename(out_html),
  output_dir    = dirname(out_html),
  knit_root_dir = getwd(),
  quiet         = FALSE)

cat(sprintf("\n[report] DONE -> %s\n", out_html))
