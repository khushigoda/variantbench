#!/usr/bin/env Rscript
# =============================================================================
# 05a_coloc.R — colocalisation of LDL fine-mapping signals with cis-QTLs
# =============================================================================
# Tests whether each fine-mapped LDL locus shares a causal variant with the
# cis-QTL of its effector gene (H4 = shared causal). Single-causal-variant
# coloc (coloc.abf), scaled to a single trait on CPU.
#
# METHOD CHOICE (Option 1 — single-signal ABF). We use coloc::coloc.abf on
# MARGINAL summary statistics from both traits over the shared cis window:
#   GWAS side : tabix-slice the locus window from the bgzipped meta stats
#               ({TRAIT}.meta.tsv.bgz) -> beta, se, eaf, N per variant.
#   QTL side  : eQTL Catalogue v2 /associations REST API, queried by the
#               specific fine-mapped molecular trait (best_mtid from the
#               coverage scan) -> beta, se, maf, an(=2N) per variant.
#   Engine    : coloc.abf(D_gwas, D_qtl), both type="quant", Wakefield ABF.
#
# WHY coloc.abf and not SuSiE-coloc (coloc.susie / coloc.bf_bf):
#   coloc.abf assumes AT MOST ONE causal variant per trait in the window. It
#   needs only marginal beta/se (no LD matrix, no per-dataset SuSiE fit, no
#   500 MB lBF download), so it runs on CPU against the light REST API. Its
#   cost is the single-causal assumption: at a locus with two independent
#   causal variants (allelic heterogeneity) it can miss a true colocalisation
#   or blur H3/H4. The multi-signal alternative (SuSiE lBFs -> coloc.bf_bf,
#   what the authors' gpu-coloc runs) relaxes that assumption but needs the
#   large per-dataset lBF matrices and a matched LD reference. That is written
#   up as a documented limitation + future refinement (Step 5 report section).
#
# TISSUE/METHOD CHOICE IS DATA-DRIVEN, NOT ASSUMED. Liver is the biologically
# obvious tissue for LDL, but GTEx liver (n=208) is the ONLY liver dataset in
# the Catalogue and fine-maps none of our effector genes. So 05_coverage_scan.py
# sweeps the WHOLE Catalogue across ALL quant methods (ge/exon/tx/txrev/
# leafcutter/microarray/aptamer) and, per (gene x quant method), names the best-
# powered dataset that fine-maps it. This script colocs each of those candidates
# so a splice-QTL signal (e.g. HMGCR via leafcutter) is tested on its own terms,
# not masked by a better-powered but non-colocalising ge dataset.
#
# INPUT  {TRAIT}.eqtl_coverage.best.tsv  (from 05_coverage_scan.py):
#   gene ensg dataset study_id study tissue quant n n_cs max_pip n_mtid best_mtid
#
# Env vars (all optional except TRAIT):
#   TRAIT LDDIR OUTDIR TARGETS
#   BEST      (coverage best.tsv; default {OUTDIR}/{TRAIT}.eqtl_coverage.best.tsv)
#   GWAS_BGZ  (tabix'd meta stats; default results/meta/{TRAIT}.meta.tsv.bgz)
#   API       (eQTL Catalogue REST base)
#   P1 P2 P12 PP4_MIN  GENES (comma list of symbols, for smoke tests)
# =============================================================================

suppressMessages({ library(data.table); library(coloc); library(jsonlite) })

env <- function(k, d="") Sys.getenv(k, unset=d)

TRAIT   <- env("TRAIT",  "LDL_C")
LDDIR   <- env("LDDIR",  "results/ld")
OUTDIR  <- env("OUTDIR", "results/coloc")
TARGETS <- env("TARGETS", sprintf("config/coloc_targets.%s.tsv", TRAIT))
BEST    <- env("BEST", file.path(OUTDIR, sprintf("%s.eqtl_coverage.best.tsv", TRAIT)))
GWAS_BGZ<- env("GWAS_BGZ", sprintf("results/meta/%s.meta.tsv.bgz", TRAIT))
API     <- env("API", "https://www.ebi.ac.uk/eqtl/api/v2")
P1      <- as.numeric(env("P1",  "1e-4"))
P2      <- as.numeric(env("P2",  "1e-4"))
P12     <- as.numeric(env("P12", "5e-6"))
PP4MIN  <- as.numeric(env("PP4_MIN", "0.8"))
GENES   <- env("GENES", "")   # optional comma-separated symbol subset (smoke tests)
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)
cat(sprintf("== 05a coloc (coloc.abf, single-signal) | trait=%s | p12=%.1e PP4>=%.2f ==\n",
            TRAIT, P12, PP4MIN))

# ---- best QTL candidate per (gene x quant method) --------------------------
if (!file.exists(BEST))
  stop("coverage best.tsv not found: ", BEST,
       "\n   run first:  TRAIT=", TRAIT, " python workflow/scripts/05_coverage_scan.py")
bt <- fread(BEST)
stopifnot(all(c("gene","ensg","dataset","study_id","tissue","quant","n","best_mtid") %in% names(bt)))
bt <- bt[dataset != "NONE" & !is.na(dataset) & best_mtid != "NA"]
if (nzchar(GENES)) bt <- bt[gene %in% strsplit(GENES, ",")[[1]]]
cat(sprintf("   %d (gene x quant) candidate(s) across %d gene(s), %d dataset(s)\n",
            nrow(bt), uniqueN(bt$ensg), uniqueN(bt$dataset)))

# ---- loci table: locus/gene -> chr, window (defines GWAS region) -----------
loci <- fread(file.path(LDDIR, sprintf("%s.loci.tsv", TRAIT)))
setnames(loci, names(loci), tolower(names(loci)))
# col 'gene' is the locus label (may be a fused name e.g. ABCG5_8); map by ensg
tgmap <- fread(TARGETS)                          # locus ensg symbol
setnames(tgmap, names(tgmap), tolower(names(tgmap)))
loci_by_locus <- merge(loci[, .(locus=gene, chr, start, end)],
                       tgmap[, .(locus, ensg)], by="locus", all.x=TRUE)
loci_by_ensg  <- loci_by_locus[!is.na(ensg)]

# ---- harmonisation key: chr:pos:sorted-alleles (orientation independent) ---
mkkey <- function(chr, pos, a1, a2) {
  chr <- sub("^chr", "", as.character(chr))
  al  <- toupper(paste(pmin(a1, a2), pmax(a1, a2), sep="_"))
  paste(chr, pos, al, sep=":")
}

# ---- GWAS side: tabix-slice the locus window from the bgzipped meta stats --
# header: SNP chr pos ea oa eaf beta se z pval neg_log10p n direction het_q het_p n_cohorts
gwas_region <- function(chr, start, end) {
  reg <- sprintf("%s:%d-%d", sub("^chr","",as.character(chr)), start, end)
  txt <- tryCatch(system2("tabix", c(shQuote(GWAS_BGZ), reg), stdout=TRUE, stderr=FALSE),
                  error=function(e) character(0))
  if (length(txt) == 0) return(data.table())
  d <- fread(text=paste(txt, collapse="\n"), header=FALSE)
  if (ncol(d) < 12) return(data.table())
  setnames(d, 1:12, c("SNP","chr","pos","ea","oa","eaf","beta","se","z","pval","nlp","n"))
  d[, key := mkkey(chr, pos, ea, oa)]
  d[is.finite(beta) & is.finite(se) & se > 0 & is.finite(eaf) & eaf > 0 & eaf < 1]
}

# ---- QTL side: eQTL Catalogue /associations for one molecular trait --------
# We query by gene_id (works for EVERY quant method) and filter to the specific
# fine-mapped molecular trait client-side. The molecular_trait_id filter cannot
# be used directly: for exon/tx/txrev/leafcutter the trait id is compound
# (e.g. ENSG....grp_1.upstream.ENST..., ENSG..._8_19955841_19956004) and the
# API's expression parser 400s on the '.'/'_' characters. gene_id returns all
# sub-features of the gene (one cis window), then we keep only best_mtid.
# Paginates by start/size; returns beta se maf N(=an/2) + orientation-free key.
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
  need <- c("beta","se","maf","an","chromosome","position","ref","alt","molecular_trait_id")
  if (!all(need %in% names(q))) return(data.table())
  # keep only the fine-mapped molecular trait (one exon/transcript/intron/gene)
  q <- q[molecular_trait_id == mtid]
  q <- q[is.finite(beta) & is.finite(se) & se > 0 &
         is.finite(maf) & maf > 0 & maf < 0.5 & is.finite(an) & an > 0]
  if (nrow(q) == 0) return(data.table())
  q[, key := mkkey(chromosome, position, ref, alt)]
  q[, N := as.integer(round(an/2))]
  q[]
}

# ---- one coloc.abf on marginal stats over the shared key set ---------------
coloc_pair <- function(gw, qt) {
  gw <- gw[!duplicated(key)]; qt <- qt[!duplicated(key)]
  shared <- intersect(gw$key, qt$key)
  if (length(shared) < 5) return(list(err=sprintf("only %d shared variants", length(shared))))
  g <- gw[match(shared, key)]; q <- qt[match(shared, key)]
  Dg <- list(beta=g$beta, varbeta=g$se^2, snp=shared, position=g$pos,
             type="quant", N=as.integer(max(g$n, na.rm=TRUE)), MAF=pmin(g$eaf, 1-g$eaf))
  Dq <- list(beta=q$beta, varbeta=q$se^2, snp=shared, position=q$position,
             type="quant", N=as.integer(max(q$N, na.rm=TRUE)), MAF=q$maf)
  res <- tryCatch(
    suppressWarnings(coloc.abf(Dg, Dq, p1=P1, p2=P2, p12=P12)),
    error=function(e) list(err=conditionMessage(e)))
  if (!is.null(res$err)) return(res)
  s <- as.data.table(as.list(res$summary))
  setnames(s, "nsnps", "nsnp", skip_absent=TRUE)
  s[, nshared := length(shared)]
  s
}

# ---------------------------------------------------------------------------
# main loop — one coloc per (gene x quant method) candidate
# ---------------------------------------------------------------------------
rows <- list()
gwas_cache <- list()   # locus -> GWAS region slice (reused across quant methods)
for (i in seq_len(nrow(bt))) {
  ensg <- bt$ensg[i]; sym <- bt$gene[i]; ds <- bt$dataset[i]
  quant <- bt$quant[i]; tissue <- bt$tissue[i]; eqtl_n <- bt$n[i]; mtid <- bt$best_mtid[i]
  lc <- loci_by_ensg[ensg == bt$ensg[i]][1]
  if (nrow(lc) == 0 || is.na(lc$chr)) {
    rows[[length(rows)+1]] <- data.table(symbol=sym, ensg, dataset=ds, quant, status="no_locus_window")
    next
  }
  locus <- lc$locus; chr <- lc$chr

  # GWAS region (cache per locus — shared across this gene's quant candidates)
  if (is.null(gwas_cache[[locus]])) gwas_cache[[locus]] <- gwas_region(chr, lc$start, lc$end)
  gw <- gwas_cache[[locus]]
  if (nrow(gw) == 0) {
    rows[[length(rows)+1]] <- data.table(locus, symbol=sym, ensg, dataset=ds, quant,
                                         status="no_gwas_region")
    next
  }

  # QTL associations for the fine-mapped molecular trait
  qt <- qtl_assoc(ds, ensg, mtid)
  if (nrow(qt) == 0) {
    rows[[length(rows)+1]] <- data.table(locus, symbol=sym, ensg, dataset=ds, quant,
                                         qtl=tissue, mtid, status="no_qtl_assoc")
    next
  }

  cc <- coloc_pair(gw, qt)
  if (!is.null(cc$err)) {
    rows[[length(rows)+1]] <- data.table(locus, symbol=sym, ensg, dataset=ds, quant,
                                         qtl=tissue, mtid, status=cc$err)
    next
  }
  cc[, `:=`(trait=TRAIT, locus=locus, symbol=sym, ensg=ensg, dataset=ds, quant=quant,
            qtl=tissue, eqtl_n=eqtl_n, mtid=mtid, status="ok")]
  rows[[length(rows)+1]] <- cc
  cat(sprintf("  %-8s %-7s %-10s %-20s n=%s  PP.H4=%.3f (H3=%.3f)  nsnp=%d\n",
              locus, sym, quant, substr(tissue,1,20), eqtl_n,
              cc$PP.H4.abf, cc$PP.H3.abf, cc$nsnp))
}

res <- rbindlist(rows, fill=TRUE)
front <- intersect(c("trait","locus","symbol","ensg","quant","qtl","dataset","mtid",
                     "eqtl_n","status","nsnp","nshared",
                     "PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf"),
                   names(res))
setcolorder(res, c(front, setdiff(names(res), front)))
out_all <- file.path(OUTDIR, sprintf("%s.coloc.tsv", TRAIT))
fwrite(res, out_all, sep="\t")

# per-gene best candidate (max PP.H4 across quant methods) + coloc call
ok <- res[status == "ok"]
if (nrow(ok)) {
  best <- ok[ok[, .I[which.max(PP.H4.abf)], by=.(ensg)]$V1]
  best[, colocalises := PP.H4.abf >= PP4MIN]
  fwrite(best, file.path(OUTDIR, sprintf("%s.coloc.best.tsv", TRAIT)), sep="\t")
  nhit <- best[colocalises == TRUE, .N]
  cat(sprintf("\n[done] %d gene(s) tested (%d gene x quant pairs) | %d colocalise (PP.H4>=%.2f)\n",
              nrow(best), nrow(ok), nhit, PP4MIN))
  cat(sprintf("       -> %s + .best.tsv\n", basename(out_all)))
  if (nhit) cat("   COLOC: ",
                paste(best[colocalises==TRUE,
                           sprintf("%s/%s(%.2f)", symbol, quant, PP.H4.abf)],
                      collapse=", "), "\n")
} else {
  cat("\n[done] no genes produced a coloc result — check statuses in", basename(out_all), "\n")
}
