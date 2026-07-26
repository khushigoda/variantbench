#!/usr/bin/env Rscript
# =====================================================================================
# 04a_ld_helpers.R — the R half of Step-4 LD prep (base R only; no packages).
#
# Split out from 04b_make_ld.sh so the R logic runs under the R conda env and the
# plink2 / download logic stays in bash under the plink env (one env per tool).
# Called by 04b with an explicit subcommand:
#
#   Rscript 04a_ld_helpers.R select    <lead.tsv> <out.loci.tsv> <W> <RAD> <curate01> <anchors.tsv> [LOCI_csv]
#   Rscript 04a_ld_helpers.R keep      <psam> <out.keep>
#   Rscript 04a_ld_helpers.R harmonize <panel.pvar> <gwas.tsv> <out_stem> <gene> <chr> <lead_pos> <W> <EDGE> <manifest>
#
# Every routine is deterministic and writes plain TSVs that plink2 (extract/LD) consumes.
# =====================================================================================
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("usage: 04a_ld_helpers.R <select|keep|harmonize> ...")
cmd  <- args[1]; a <- args[-1]

# ------------------------------------------------------------------ subcommand: select ---
# Eligibility screen (both modes) + optional CAD-effector curation. Writes the loci table
# whose chromosome column drives which per-chr panel files 04b downloads.
if (cmd == "select") {
  lead <- read.delim(a[1], stringsAsFactors = FALSE)
  out  <- a[2]; W <- as.numeric(a[3]); RAD <- as.numeric(a[4])
  curate <- a[5] == "1"; anchors_f <- a[6]
  sel_env <- if (length(a) >= 7) a[7] else ""     # optional comma-list (sample test)
  n0 <- nrow(lead)
  # --- eligibility screen (applies in BOTH modes; log attrition) ---
  lead <- lead[lead$maf_class != "unknown", ]                                   ; n1 <- nrow(lead)
  lead <- lead[!(lead$chr == 6 & lead$pos >= 25e6 & lead$pos <= 34e6), ]        ; n2 <- nrow(lead)
  lead <- lead[lead$maf_class != "rare", ]                                      ; n3 <- nrow(lead)
  cat(sprintf("    leads: %d -> drop unknown-freq -> %d -> drop MHC -> %d -> drop rare -> %d eligible\n",
              n0, n1, n2, n3))
  rows <- list()
  if (!curate) {
    # CURATE=0: every eligible lead is its own locus (whole lead list / other metabolites)
    for (i in seq_len(nrow(lead))) rows[[length(rows)+1]] <-
      data.frame(gene = paste0("lead_", lead$SNP[i]), chr = lead$chr[i],
                 lead_snp = lead$SNP[i], lead_pos = lead$pos[i],
                 ea = lead$ea[i], oa = lead$oa[i], z = lead$z[i],
                 nlp = lead$neg_log10p[i], stringsAsFactors = FALSE)
  } else {
    # CURATE=1: strongest eligible lead near each gene anchor (anchors from external TSV)
    genes <- read.table(anchors_f, header = FALSE, stringsAsFactors = FALSE,
                        col.names = c("gene","chr","anchor"), comment.char = "#")
    want <- if (nzchar(sel_env)) strsplit(sel_env, ",")[[1]] else genes$gene
    for (g in want) {
      gi <- genes[genes$gene == g, ]
      if (!nrow(gi)) { cat(sprintf("    [skip] unknown gene '%s'\n", g)); next }
      cand <- lead[lead$chr == gi$chr & abs(lead$pos - gi$anchor) <= RAD, ]
      if (!nrow(cand)) { cat(sprintf("    [skip] %s: no eligible lead within %.0fkb of anchor\n",
                                      g, RAD/1e3)); next }
      cand <- cand[order(-cand$neg_log10p), ][1, ]   # strongest eligible lead near anchor
      rows[[length(rows)+1]] <-
        data.frame(gene = g, chr = gi$chr, lead_snp = cand$SNP, lead_pos = cand$pos,
                   ea = cand$ea, oa = cand$oa, z = cand$z, nlp = cand$neg_log10p,
                   stringsAsFactors = FALSE)
    }
  }
  loci <- do.call(rbind, rows)
  loci$start <- pmax(1, loci$lead_pos - W)
  loci$end   <- loci$lead_pos + W
  write.table(loci, out, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("    selected %d loci -> %s\n", nrow(loci), out))
  print(loci[, c("gene","chr","lead_snp","lead_pos","nlp")], row.names = FALSE)

# -------------------------------------------------------------------- subcommand: keep ---
# EUR keep file from .psam. Schema varies (#FID+IID vs #IID; SuperPop col name) — detect it.
} else if (cmd == "keep") {
  psam_f <- a[1]; keep <- a[2]
  psam <- readLines(psam_f)
  hdr_i <- grep("^#", psam); hdr_i <- if(length(hdr_i)) hdr_i[1] else 1
  hdr <- strsplit(sub("^#", "", psam[hdr_i]), "\t")[[1]]
  dat <- read.table(text = psam[(hdr_i+1):length(psam)], sep = "\t",
                    stringsAsFactors = FALSE, col.names = hdr, check.names = FALSE)
  sp <- grep("^superpop$|^super_pop$|^pop$|population", hdr, ignore.case = TRUE, value = TRUE)
  spcol <- NA
  for (c in sp) if (any(dat[[c]] == "EUR")) { spcol <- c; break }
  if (is.na(spcol)) stop("no SuperPop column with value 'EUR' found in .psam")
  eur <- dat[dat[[spcol]] == "EUR", , drop = FALSE]
  idcols <- intersect(c("FID","IID"), hdr); if(!length(idcols)) idcols <- "IID"
  out <- eur[, idcols, drop = FALSE]
  writeLines(paste0("#", paste(idcols, collapse = "\t")), keep)
  write.table(out, keep, sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = FALSE, append = TRUE)
  cat(sprintf("    EUR samples: %d   (id cols: %s;  superpop col: %s)\n",
              nrow(eur), paste(idcols, collapse="+"), spcol))

# --------------------------------------------------------------- subcommand: harmonize ---
# Match panel<->GWAS on chr:pos + allele pair; allow ref/alt swap (flip z sign); drop
# strand-ambiguous palindromes. Writes .extract.ids + .z.tsv (position-ordered) + manifest row.
} else if (cmd == "harmonize") {
  pvar_f<-a[1]; gwas_f<-a[2]; stem<-a[3]; gene<-a[4]; chr<-a[5]
  lead_pos<-as.numeric(a[6]); W<-as.numeric(a[7]); EDGE<-as.numeric(a[8]); manifest<-a[9]
  pv <- read.table(pvar_f, sep="\t", stringsAsFactors=FALSE, comment.char="")
  names(pv)[1:5] <- c("CHROM","POS","ID","REF","ALT")   # pvar: #CHROM POS ID REF ALT (+maybe more)
  gw <- read.delim(gwas_f, stringsAsFactors=FALSE)
  n_panel <- nrow(pv); n_gwas <- nrow(gw)
  key <- function(chr,pos) paste(chr,pos,sep=":")
  gw$k <- key(gw$chr, gw$pos); pv$k <- key(pv$CHROM, pv$POS)
  m <- merge(pv, gw, by="k")                             # match on chr:pos
  n_match_pos <- nrow(m)
  palin <- function(a,b){ p<-paste0(toupper(a),toupper(b)); p %in% c("AT","TA","CG","GC") }
  aligned_z <- rep(NA_real_, nrow(m)); keepflag <- logical(nrow(m)); swapped <- 0L; ndrop_pal <- 0L
  for (i in seq_len(nrow(m))) {
    REF<-toupper(m$REF[i]); ALT<-toupper(m$ALT[i]); EA<-toupper(m$ea[i]); OA<-toupper(m$oa[i])
    if (palin(REF,ALT)) { ndrop_pal<-ndrop_pal+1L; next }   # strand-ambiguous -> drop
    if (EA==ALT && OA==REF)      { aligned_z[i]<- m$z[i];  keepflag[i]<-TRUE }  # GWAS EA==ALT: z as-is
    else if (EA==REF && OA==ALT) { aligned_z[i]<- -m$z[i]; keepflag[i]<-TRUE; swapped<-swapped+1L }  # swap -> flip
    # else: allele pair mismatch -> drop (keepflag stays FALSE)
  }
  m$aligned_z <- aligned_z; keep <- m[keepflag, ]
  keep <- keep[order(as.numeric(keep$POS)), ]            # POSITION order == pvar/LD order
  n_final <- nrow(keep)
  edge <- "no"
  if (n_final > 0) {
    top <- keep[which.max(keep$nlp), ]
    if (abs(top$POS - (lead_pos - W)) <= EDGE || abs(top$POS - (lead_pos + W)) <= EDGE) edge <- "YES"
  }
  writeLines(keep$ID, paste0(stem, ".extract.ids"))       # for plink2 --extract
  write.table(data.frame(ID=keep$ID, POS=keep$POS, z=keep$aligned_z),
              paste0(stem, ".z.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  status <- if (n_final >= 2) "OK" else "TOO_FEW_VARIANTS"
  cat(sprintf("%s\t%s\t%s\t%.0f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n",
      gene, chr, "", lead_pos, 2*W, n_panel, n_gwas, n_match_pos, swapped, ndrop_pal, n_final, edge, status),
      file=manifest, append=TRUE)
  cat(sprintf("       panel=%d gwas=%d pos-match=%d swapped=%d palindrome-drop=%d final=%d edge=%s\n",
      n_panel, n_gwas, n_match_pos, swapped, ndrop_pal, n_final, edge))

} else stop(sprintf("unknown subcommand '%s' (select|keep|harmonize)", cmd))
