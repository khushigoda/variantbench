#!/usr/bin/env Rscript
# =====================================================================
# 04e — Fine-mapping visualization (report figures)
# =====================================================================
# Renders the figures for the take-home report from the persisted 04c/04d
# tables ONLY — never re-fits SuSiE, never re-queries VEP. Reads:
#   results/finemap/{trait}.finemap.tsv     per-variant POS/z/pip/cs95/r2_lead
#   results/finemap/{trait}.cs_summary.tsv  per-CS pass_filter + purity (95%)
#   results/finemap/{trait}.annot.tsv       consequence per CS variant (04d)
#
# Produces three deliverables:
#   1. per-locus regional plots  plots/{trait}.{gene}.regional.png
#        two stacked tracks sharing the x-axis (LocusZoom idiom):
#        top = -log10(P) coloured by r2-to-lead; bottom = PIP with CS variants
#        highlighted and MISSENSE variants labelled (gene + aa_change) — the
#        figure that visually answers "which fine-mapped variant is missense".
#   2. cross-locus QC summary    {trait}.finemap_summary.png
#        per-locus CS counts (raw vs pass-filter) + convergence + kriging_s.
#   3. consequence + missense    {trait}.consequence_summary.png
#        stacked consequence breakdown (coding-vs-regulatory headline) and the
#        missense hit list with SIFT/PolyPhen.
#
# House style mirrors 03_leadvars.R: theme_bw(11), geom_text(check_overlap),
# ggsave(w,h,dpi) — ggplot2 + data.table only, no extra packages.
#
# Env (mirrors 04c/04d): TRAIT FMDIR LDDIR OUTDIR LOCI(optional subset)
suppressMessages({library(data.table); library(ggplot2)})

TRAIT  <- Sys.getenv("TRAIT",  "LDL_C")
FMDIR  <- Sys.getenv("FMDIR",  "results/finemap")
LDDIR  <- Sys.getenv("LDDIR",  "results/ld")
OUTDIR <- Sys.getenv("OUTDIR", "results/finemap")
LOCI   <- Sys.getenv("LOCI",   "")               # optional CSV gene subset
PLOTDIR <- file.path(OUTDIR, "plots")
dir.create(PLOTDIR, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("== 04e fine-mapping plots | trait=%s ==\n", TRAIT))

fm  <- fread(file.path(FMDIR, paste0(TRAIT, ".finemap.tsv")))
css <- fread(file.path(FMDIR, paste0(TRAIT, ".cs_summary.tsv")))
apath <- file.path(FMDIR, paste0(TRAIT, ".annot.tsv"))
ann <- if (file.exists(apath)) fread(apath) else data.table()

# -log10(P) from z (two-sided); robust to huge |z| (values that overflow -> cap)
fm[, nlp := pmin(-(pnorm(abs(z), lower.tail = FALSE, log.p = TRUE) + log(2)) / log(10), 700)]

genes <- if (nzchar(LOCI)) strsplit(LOCI, ",")[[1]] else sort(unique(fm$gene))

# ---------------------------------------------------------------- 1. regional
band_lvls <- c("0.8-1.0","0.6-0.8","0.4-0.6","0.2-0.4","0.0-0.2","lead")
band_cols <- c("lead"="#7B3294","0.8-1.0"="#D7191C","0.6-0.8"="#FDAE61",
               "0.4-0.6"="#FFFFBF","0.2-0.4"="#ABD9E9","0.0-0.2"="#2C7BB6")
r2band <- function(r2, is_lead) {
  b <- cut(r2, c(-Inf,0.2,0.4,0.6,0.8,Inf),
           labels=c("0.0-0.2","0.2-0.4","0.4-0.6","0.6-0.8","0.8-1.0"))
  b <- as.character(b); b[is.na(b)] <- "0.0-0.2"; b[is_lead==1] <- "lead"; b
}

regional_one <- function(g) {
  d <- fm[gene == g]
  if (!nrow(d)) return(invisible(NULL))
  d[, band := factor(r2band(r2_lead, is_lead), levels = band_lvls)]
  d[, mb := POS/1e6]
  d[, incs := cs95 != 0 & !is.na(cs95)]
  # missense variants at this locus (from 04d), for callout labels
  ms <- if (nrow(ann)) ann[gene == g & is_missense == TRUE] else data.table()

  top <- ggplot(d, aes(mb, nlp, colour = band)) +
    geom_point(size = 0.7, alpha = 0.8) +
    scale_colour_manual(values = band_cols, name = expression(r^2~"to lead"),
                        drop = FALSE) +
    labs(x = NULL, y = expression(-log[10](italic(P))),
         title = sprintf("%s — %s locus", TRAIT, g)) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "right",
          plot.title = element_text(face = "bold"))

  bot <- ggplot(d, aes(mb, pip)) +
    geom_point(aes(colour = incs), size = 0.7, alpha = 0.7) +
    scale_colour_manual(values = c("FALSE"="grey70","TRUE"="#D7191C"),
                        labels = c("not in CS","95% CS"), name = NULL) +
    scale_y_continuous(limits = c(0, 1.08), breaks = c(0,0.5,1)) +
    labs(x = sprintf("Chr %s position (Mb)", d$chr[1]), y = "PIP") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  if (nrow(ms)) {
    ms[, mb := POS/1e6]
    bot <- bot +
      geom_point(data = ms, aes(mb, pip), colour = "#1A9641", size = 2.2,
                 shape = 18, inherit.aes = FALSE) +
      geom_text(data = ms, aes(mb, pip,
                 label = paste0(gene_symbol, " ", aa_change)),
                 inherit.aes = FALSE, size = 2.9, fontface = "italic",
                 vjust = -0.8, hjust = 0.5, colour = "#1A9641",
                 check_overlap = TRUE)
  }
  # stack the two tracks without patchwork: render to one column via a shared
  # facet on a long-melt is messy for differing y-scales, so save each track
  # and combine with grid layout via a 2-row viewport.
  png(file.path(PLOTDIR, sprintf("%s.%s.regional.png", TRAIT, g)),
      width = 2000, height = 1500, res = 200)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 1, heights = c(1.15, 1))))
  print(top, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(bot, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  grid::popViewport()
  dev.off()
  invisible(TRUE)
}
for (g in genes) regional_one(g)
cat(sprintf("   regional plots: %d loci -> %s/\n", length(genes), PLOTDIR))

# ---------------------------------------------------------------- 2. QC summary
qc <- fread(file.path(FMDIR, paste0(TRAIT, ".susie_qc.tsv")))
cs95 <- css[abs(coverage - 0.95) < 1e-9]
cnt <- cs95[, .(n_raw = .N, n_pass = sum(pass_filter == TRUE)), by = gene]
cnt <- merge(cnt, qc[, .(gene, converged, kriging_s)], by = "gene", all.x = TRUE)
cnt <- cnt[order(-n_pass)]
cnt[, gene := factor(gene, levels = gene)]
m <- melt(cnt, id.vars = c("gene","converged","kriging_s"),
          measure.vars = c("n_raw","n_pass"),
          variable.name = "kind", value.name = "n")
m[, kind := factor(kind, levels=c("n_raw","n_pass"),
                   labels=c("all CS","pass filter (genome-wide sig)"))]

p_qc <- ggplot(m, aes(gene, n, fill = kind)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(data = cnt, aes(gene, y = -0.4,
            label = ifelse(converged=="TRUE" | converged==TRUE, "conv", "flag")),
            inherit.aes = FALSE, size = 2.6,
            colour = ifelse(cnt$converged=="TRUE" | cnt$converged==TRUE, "grey40","#D7191C")) +
  scale_fill_manual(values = c("all CS"="#ABD9E9","pass filter (genome-wide sig)"="#2C7BB6"),
                    name = NULL) +
  labs(x = NULL, y = "credible sets (95%)",
       title = sprintf("%s — fine-mapping QC across %d loci", TRAIT, nrow(cnt)),
       subtitle = "bars = CS count; 'flag' = did not formally converge (out-of-sample LD)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        legend.position = "top", plot.title = element_text(face="bold"),
        panel.grid.major.x = element_blank())
ggsave(file.path(OUTDIR, paste0(TRAIT, ".finemap_summary.png")), p_qc,
       width = 9, height = 5, dpi = 150)
cat("   QC summary -> ", file.path(OUTDIR, paste0(TRAIT, ".finemap_summary.png")), "\n")

# ---------------------------------------------------------------- 3. consequence + missense
if (nrow(ann)) {
  # consequence breakdown (coding-vs-regulatory headline)
  coding <- c("missense_variant","synonymous_variant","stop_gained","stop_lost",
              "start_lost","splice_donor_variant","splice_acceptor_variant",
              "splice_region_variant","splice_polypyrimidine_tract_variant")
  ann[, class := ifelse(consequence %in% coding, "coding", "regulatory / non-coding")]
  cb <- ann[, .N, by = .(consequence, class)][order(class, -N)]
  cb[, consequence := factor(consequence, levels = rev(cb$consequence))]
  p_cons <- ggplot(cb, aes(consequence, N, fill = class)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = N), hjust = -0.2, size = 3) +
    coord_flip() +
    scale_fill_manual(values = c("coding"="#D7191C","regulatory / non-coding"="#2C7BB6"),
                      name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = NULL, y = "95% CS variants",
         title = sprintf("%s — consequence of fine-mapped variants", TRAIT),
         subtitle = sprintf("%d CS variants; only %d missense — most causal signal is regulatory",
                            nrow(ann), sum(ann$is_missense == TRUE))) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top", plot.title = element_text(face="bold"),
          panel.grid.major.y = element_blank())
  ggsave(file.path(OUTDIR, paste0(TRAIT, ".consequence_summary.png")), p_cons,
         width = 8, height = 5, dpi = 150)
  cat("   consequence summary -> ",
      file.path(OUTDIR, paste0(TRAIT, ".consequence_summary.png")), "\n")

  # missense hit list as a table-figure
  ms <- ann[is_missense == TRUE][order(gene, -pip)]
  if (nrow(ms)) {
    ms[, lab_sift := ifelse(sift_pred=="", "-", sprintf("%s (%.2g)", sift_pred, sift_score))]
    ms[, lab_pp   := ifelse(polyphen_pred=="", "-", sprintf("%s (%.2g)", polyphen_pred, polyphen_score))]
    tab <- ms[, .(Gene=gene_symbol, Variant=ID, Change=aa_change,
                  PIP=sprintf("%.2f", pip), CS=cs,
                  SIFT=lab_sift, PolyPhen=lab_pp,
                  GT=ifelse(is_known_causal==TRUE,"*",""))]
    # render as a simple grid table
    ncol_t <- ncol(tab); nrow_t <- nrow(tab)
    gt <- ggplot() + theme_void() +
      xlim(0.5, ncol_t + 0.5) + ylim(0.5, nrow_t + 1.5) +
      labs(title = sprintf("%s — missense variants shown by fine-mapping", TRAIT),
           subtitle = "* = known drug-target positive control (PCSK9 R46L)") +
      theme(plot.title = element_text(face="bold", size=12),
            plot.subtitle = element_text(size=9))
    hdr <- names(tab)
    for (j in seq_len(ncol_t))
      gt <- gt + annotate("text", x=j, y=nrow_t+1, label=hdr[j],
                          fontface="bold", size=3.2, hjust=0.5)
    for (i in seq_len(nrow_t)) for (j in seq_len(ncol_t))
      gt <- gt + annotate("text", x=j, y=nrow_t+1-i,
                          label=as.character(tab[[j]][i]), size=3, hjust=0.5)
    ggsave(file.path(OUTDIR, paste0(TRAIT, ".missense_table.png")), gt,
           width = 11, height = 1 + 0.5*nrow_t, dpi = 150)
    cat("   missense table -> ",
        file.path(OUTDIR, paste0(TRAIT, ".missense_table.png")), "\n")
  }
} else {
  cat("   (annot.tsv absent — skipping consequence/missense figures)\n")
}

cat(sprintf("[done] %d regional + summary figures written under %s\n",
            length(genes), OUTDIR))
