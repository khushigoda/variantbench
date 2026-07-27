# Step 5 (colocalization) — design notes (in progress)

**Status: design under active revision — not yet run.** This document records the method
reasoning so the choices are auditable; the final architecture is being settled with the
hiring team's brief in mind ("we care more about method choices and reasoning than results").

## Goal
For each of the 11 LDL fine-mapped effector loci, test whether the LDL GWAS signal and a
**molecular-QTL** signal share the same causal variant — i.e. is the LDL association mediated
through a measurable molecular phenotype (expression, splicing, protein level) of the effector
gene? This is trait–QTL colocalization against the **eQTL Catalogue**.

Target genes (`config/coloc_targets.LDL_C.tsv`): PCSK9, APOB, ABCG5, ABCG8, HMGCR, NPC1L1,
LPL, ABCA1, APOA5, CETP, LDLR, APOE (12 genes across 11 loci; ABCG5+ABCG8 share one locus).

## What "molecular QTL" scopes in
A molecular QTL associates a variant with an **intermediate molecular phenotype**, not the
organismal trait. The eQTL Catalogue quantifies these in several ways, and the scope is *all*
of them, not just total gene expression:
- `ge` — gene expression (classic eQTL)
- `exon` — per-exon expression
- `tx` — transcript-level (Salmon) abundance
- `txrev` — transcript-usage / splicing events
- `leafcutter` — intron-excision splicing QTLs (sQTL)
- `microarray` — array-based expression
- protein QTLs (aptamer/Somalogic) where present

Restricting to `ge` alone would miss, e.g., **HMGCR**, where the LDL signal acts through a
splice-altering variant (exon skipping) — visible to `leafcutter`/`txrev` but not to a
total-expression scan. Trait–trait coloc against the CAD GWAS is a separate, optional bonus,
outside the molecular-QTL definition.

## Two candidate architectures (the open decision)
The take-home authors run tissue-agnostic discovery at biobank scale on GPU (a reimplementation
of coloc's Bayes-factor engine). This project scales that down to a hypothesis-driven, targeted
coloc on CPU using the **same underlying Bayes-factor model** — a deliberate, stated scaling-down.
Two faithful ways to do it, differing in single- vs multi-signal resolution:

1. **Region `/associations` API + `coloc.abf` (single-signal), all quant methods.** Pull
   region-filtered summary statistics (β, se, MAF, N) per gene per dataset via the Catalogue's
   paged API — the maintainer-recommended access pattern (their API v2 tutorial). Run `coloc.abf`,
   which assumes ≤1 causal variant per region. Scales cleanly across all quant methods and the
   whole catalogue; the cost is the single-causal-variant assumption, which can under-call at
   loci with allelic heterogeneity (CETP, APOE, LDLR).
2. **Option 1 + a SuSiE `coloc.bf_bf` multi-signal cross-check** on the genes that have a
   relevant credible set, reusing our Step-4 SuSiE fits on the GWAS side and the Catalogue's
   SuSiE lBF files on the QTL side. Multi-signal resolves several independent signals per locus
   and colocalizes them pairwise — closest to the authors' method, at the cost of the heavier
   FTP lBF files and CS-coverage limits.

**Data feasibility confirmed for both.** GWAS side: `results/meta/LDL_C.meta.tsv.bgz` carries
β/se/eaf/N per variant and is tabix-indexed for region queries. QTL side: the `/associations`
API returns β/se/maf/an; the FTP SuSiE lBF files supply the multi-signal inputs.

## Single-signal vs multi-signal, briefly
Both test the same H4 hypothesis (shared causal variant). `coloc.abf` assumes one causal variant
per region and compares each trait's single posterior peak. SuSiE coloc first fine-maps each
side into independent credible sets, then colocalizes every GWAS signal against every QTL signal
pairwise — so at a locus with two independent causal variants it can tell which one is shared.
For the lipid loci with known allelic heterogeneity this distinction is substantive.

## Key power finding (why the naive choice fails)
The biologically-obvious "just use liver" choice returns **zero** hits: GTEx liver (QTD000266,
n=208) is the only liver expression dataset in the Catalogue and fine-maps **none** of the 12
effector genes at that sample size — a real eQTL-power limitation (eQTL discovery power ≪ GWAS
power), not a bug. This is itself a reportable result and forces the data-driven, pan-dataset
tissue choice rather than a single guessed tissue.

## Scripts (staged, pending architecture lock)
- `05_coverage_scan.py` — reads the small credible-set files across gene-level datasets to pick,
  per gene, the best-powered dataset that fine-maps it (feeds the FTP/SuSiE architecture).
- `05a_coloc.R` — per-gene multi-dataset coloc engine.
These will be finalized/replaced once the region-API vs FTP architecture is locked with the team.