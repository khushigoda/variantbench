#!/usr/bin/env Rscript
# =============================================================================
# Step 5b — colocalization figures (regional stacked plots + PP overview)
# =============================================================================
# Recreates the reference paper's Supp Fig 7 style regional colocalization plot
# (GWAS panel stacked over the molecular-QTL panel, shared genomic x-axis, both
# leads marked, PP.H4 annotated) for every locus that colocalizes, plus a
# single overview panel decomposing PP.H4 (shared causal variant) vs PP.H3
# (distinct causal variants) across all tested effector genes.
#
# We do NOT copy the paper's 3-panel figure: that study has a disease (CAD) GWAS
# layer on top; this replication is trait-vs-molecular-QTL only, so the honest
# analogue is 2 panels (LDL-C GWAS over molecular QTL).
#
# The coloc engine (05a_coloc.R) persists only the summary tables, not the
# per-variant regional data — so we re-slice the GWAS window from the tabix'd
# meta stats and re-fetch the QTL associations from the eQTL Catalogue REST API,
# reusing the IDENTICAL request pattern as 05a_coloc.R.
#
# Env vars (all optional except TRAIT):
#   TRAIT LDDIR OUTDIR TARGETS
#   BESTCOLOC (coloc.best.tsv; default {OUTDIR}/{TRAIT}.coloc.best.tsv)
#   GWAS_BGZ  (tabix'd meta stats; default results/meta/{TRAIT}.meta.tsv.bgz)
#   API       (eQTL Catalogue REST base)
#   PP4_MIN   (colocalization call threshold; default 0.8)
# Outputs:
#   {OUTDIR}/plots/{TRAIT}.{gene}.regional.png   (one per colocalizing gene)
#   {OUTDIR}/{TRAIT}.coloc_overview.png          (PP.H4 vs PP.H3 across genes)
# =============================================================================

suppressMessages({ library(data.table); library(jsonlite)
                   library(ggplot2); library(cowplot) })

env <- function(k, d="") Sys.getenv(k, unset=d)

TRAIT    <- env("TRAIT",  "LDL_C")
LDDIR    <- env("LDDIR",  "results/ld")
OUTDIR   <- env("OUTDIR", "results/coloc")
TARGETS  <- env("TARGETS", sprintf("config/coloc_targets.%s.tsv", TRAIT))
BESTCOLOC<- env("BESTCOLOC", file.path(OUTDIR, sprintf("%s.coloc.best.tsv", TRAIT)))
GWAS_BGZ <- env("GWAS_BGZ", sprintf("results/meta/%s.meta.tsv.bgz", TRAIT))
API      <- env("API", "https://www.ebi.ac.uk/eqtl/api/v2")
PP4MIN   <- as.numeric(env("PP4_MIN", "0.8"))
PLOTDIR  <- file.path(OUTDIR, "plots")
dir.create(PLOTDIR, recursive=TRUE, showWarnings=FALSE)
cat(sprintf("== 05b coloc plots | trait=%s | PP4>=%.2f ==\n", TRAIT, PP4MIN))

# ---- harmonisation key: chr:pos:sorted-alleles (orientation independent) ---
mkkey <- function(chr, pos, a1, a2) {
  chr <- sub("^chr", "", as.character(chr))
  al  <- toupper(paste(pmin(a1, a2), pmax(a1, a2), sep="_"))
  paste(chr, pos, al, sep=":")
}

# ---- GWAS side: tabix-slice the locus window (identical to 05a_coloc.R) -----
gwas_region <- function(chr, start, end) {
  reg <- sprintf("%s:%d-%d", sub("^chr","",as.character(chr)), start, end)
  txt <- tryCatch(system2("tabix", c(shQuote(GWAS_BGZ), reg), stdout=TRUE, stderr=FALSE),
                  error=function(e) character(0))
  if (length(txt) == 0) return(data.table())
  d <- fread(text=paste(txt, collapse="\n"), header=FALSE)
  if (ncol(d) < 12) return(data.table())
  setnames(d, 1:12, c("SNP","chr","pos","ea","oa","eaf","beta","se","z","pval","nlp","n"))
  d[is.finite(pos) & is.finite(nlp)]
}

# ---- QTL side: eQTL Catalogue /associations (identical to 05a_coloc.R) ------
qtl_assoc <- function(dataset, ensg, mtid) {
  out <- list(); start <- 0L; size <- 1000L
  repeat {
    url <- sprintf("%s/datasets/%s/associations?gene_id=%s&size=%d&start=%d",
                   API, dataset, ensg, size, start)
    d <- tryCatch(fromJSON(url), error=function(e) NULL)
    if (is.null(d) || length(d) == 0 || nrow(d) == 0) break
    out[[length(out)+1]] <- as.data.table(d)
    if (nrow(d) < size) break
    start <- start + size
  }
  if (!length(out)) return(data.table())
  q <- rbindlist(out, fill=TRUE)
  need <- c("beta","se","position","pvalue","molecular_trait_id")
  if (!all(need %in% names(q))) return(data.table())
  q <- q[molecular_trait_id == mtid]
  for (c in c("beta","se","position","pvalue")) q[[c]] <- as.numeric(q[[c]])
  q <- q[is.finite(position) & is.finite(pvalue) & pvalue > 0]
  if (nrow(q) == 0) return(data.table())
  q[, nlp := -log10(pmax(pvalue, 1e-300))]
  q[]
}

# ---- inputs ----------------------------------------------------------------
best <- fread(BESTCOLOC)
loci <- fread(file.path(LDDIR, sprintf("%s.loci.tsv", TRAIT)))
setnames(loci, names(loci), tolower(names(loci)))
tg   <- fread(TARGETS); setnames(tg, names(tg), tolower(names(tg)))
# locus window per gene symbol (loci$gene may be a fused locus label e.g. ABCG5_8)
loc_by_ensg <- merge(loci[, .(locus=gene, chr, start, end)],
                     tg[, .(locus, ensg)], by="locus", all.x=TRUE)

pp4col <- grep("PP.H4", names(best), fixed=TRUE, value=TRUE)[1]
pp3col <- grep("PP.H3", names(best), fixed=TRUE, value=TRUE)[1]
best[, coloc := get(pp4col) >= PP4MIN]

# =============================================================================
# 1. Regional stacked plots — one per colocalizing gene
# =============================================================================
regional_one <- function(sym) {
  r <- best[symbol == sym][1]
  w <- loc_by_ensg[ensg == r$ensg][1]
  if (is.na(w$chr)) return(NULL)
  g <- gwas_region(w$chr, w$start, w$end)
  q <- qtl_assoc(r$dataset, r$ensg, r$mtid)
  if (!nrow(g) || !nrow(q)) return(NULL)
  q <- q[position >= w$start & position <= w$end]      # clip QTL to GWAS window
  gl <- g[which.max(nlp)]; ql <- q[which.max(nlp)]
  xr <- range(g$pos)/1e6

  base <- theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(),
          plot.margin = margin(2, 6, 2, 6))

  pG <- ggplot(g, aes(pos/1e6, nlp)) +
    geom_point(size = 0.9, colour = "#9ecae1", alpha = 0.8) +
    geom_point(data = gl, size = 2.4, colour = "#08519c") +
    geom_text(data = gl, aes(label = SNP), vjust = -0.8, size = 2.6, colour = "#08519c") +
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed",
               colour = "grey60", linewidth = 0.3) +
    scale_x_continuous(limits = xr) +
    labs(y = expression(-log[10]*italic(P)),
         title = "LDL-C GWAS") + base +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
          plot.title = element_text(size = 8, face = "plain"))

  pQ <- ggplot(q, aes(position/1e6, nlp)) +
    geom_point(size = 0.9, colour = "#fdae6b", alpha = 0.85) +
    geom_point(data = ql, size = 2.4, colour = "#d94801") +
    scale_x_continuous(limits = xr) +
    labs(x = sprintf("chromosome %s position (Mb)", sub("^chr","",as.character(w$chr))),
         y = expression(-log[10]*italic(P)),
         title = sprintf("%s %s-QTL (%s)", sym, r$quant, r$qtl)) + base +
    theme(plot.title = element_text(size = 8, face = "plain"))

  # PP.H4 annotation banner between titles
  ttl <- ggdraw() + draw_label(
    sprintf("%s   PP.H4 = %.3f", sym, r[[pp4col]]),
    fontface = "italic", size = 10, x = 0.02, hjust = 0)
  g_stack <- plot_grid(ttl, pG, pQ, ncol = 1, rel_heights = c(0.10, 1, 1.05),
                       align = "v", axis = "lr")
  fn <- file.path(PLOTDIR, sprintf("%s.%s.regional.png", TRAIT, sym))
  ggsave(fn, g_stack, width = 5.6, height = 4.6, dpi = 150,
         device = grDevices::png, type = "cairo", bg = "white")
  # Convert to WebP (60-70% smaller)
  if (nzchar(Sys.which("cwebp"))) {
    webp_fn <- sub("\\.png$", ".webp", fn)
    system(sprintf("cwebp -q 85 '%s' -o '%s' 2>/dev/null && rm '%s'", fn, webp_fn, fn))
    fn <- webp_fn
  }
  cat(sprintf("   %-8s GWAS-lead %s@%d  QTL-peak@%d  d=%.1fkb  -> %s\n",
              sym, gl$SNP, gl$pos, ql$position, abs(gl$pos-ql$position)/1e3, basename(fn)))
  fn
}

coloc_genes <- best[coloc == TRUE][order(-get(pp4col))]$symbol
cat(sprintf("-- %d colocalizing gene(s): %s\n",
            length(coloc_genes), paste(coloc_genes, collapse=", ")))
for (sym in coloc_genes) try(regional_one(sym), silent=FALSE)

# =============================================================================
# 2. Overview — PP.H4 (shared signal) vs PP.H3 (distinct signal) per gene
# =============================================================================
ov <- best[order(get(pp4col))]
ov[, label := sprintf("%s (%s/%s)", symbol, quant, qtl)]
ov[, label := factor(label, levels = label)]
m <- melt(ov[, .(label, H4 = get(pp4col), H3 = get(pp3col))],
          id.vars = "label", variable.name = "hyp", value.name = "pp")
m[, hyp := factor(hyp, levels = c("H4","H3"),
                  labels = c("PP.H4 (shared causal variant)",
                             "PP.H3 (distinct causal variants)"))]

pOv <- ggplot(m, aes(pp, label, colour = hyp)) +
  geom_line(aes(group = label), colour = "grey80", linewidth = 0.5) +
  geom_point(size = 2.6) +
  geom_vline(xintercept = PP4MIN, linetype = "dashed",
             colour = "grey50", linewidth = 0.3) +
  scale_colour_manual(values = c("PP.H4 (shared causal variant)"    = "#c0392b",
                                 "PP.H3 (distinct causal variants)" = "#7f8c8d"),
                      name = NULL) +
  scale_x_continuous(limits = c(0, 1.02), breaks = c(0, 0.25, 0.5, 0.8, 1)) +
  labs(x = "Posterior probability",
       y = NULL,
       title = "LDL-C effector loci: colocalization with molecular QTLs") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(size = 10, face = "plain"))

ov_fn <- file.path(OUTDIR, sprintf("%s.coloc_overview.png", TRAIT))
ggsave(ov_fn, pOv, width = 6.4, height = 4.6, dpi = 150,
       device = grDevices::png, type = "cairo", bg = "white")
# Convert to WebP (60-70% smaller)
if (nzchar(Sys.which("cwebp"))) {
  webp_fn <- sub("\\.png$", ".webp", ov_fn)
  system(sprintf("cwebp -q 85 '%s' -o '%s' 2>/dev/null && rm '%s'", ov_fn, webp_fn, ov_fn))
  ov_fn <- webp_fn
}
cat(sprintf("-- overview -> %s (%d genes, %d coloc)\n",
            basename(ov_fn), nrow(best), sum(best$coloc)))
cat("[done] 05b coloc plots\n")
