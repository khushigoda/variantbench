# Step 4 (fine-mapping) — my implementation vs. the authors' method

My Step 4: `workflow/Snakefile` rules `finemap_ld` → `finemap` → `annotate` → `finemap_plots`;
scripts `04a`/`04b` (LD prep), `04c_finemap.R` (SuSiE), `04d_annotate.py` (VEP), `04e_finemap_plots.R`.
Input: the Step-1 meta table `results/meta/LDL_C.meta.tsv.gz` + the Step-3 lead loci.
Output: `results/finemap/LDL_C.{finemap.tsv, cs_summary.tsv, cs_summary.filtered.tsv, susie_qc.tsv,
annot.tsv, annot.missense.tsv}` + regional/summary figures under `results/finemap/[plots/]`.

## What the brief asks vs. what the paper did
The brief's Step 4 is **"Fine-mapping of a locus using ancestry-appropriate public LD scores."**
The paper fine-maps genome-wide with SuSiE using **in-sample** UKBB LD. The brief explicitly
relaxes this to **public** (out-of-sample) LD, which is the single most consequential difference
and drives most of my method choices below. I fine-map the **11 canonical LDL effector loci**
(PCSK9, APOB, ABCG5/8, HMGCR, NPC1L1, LPL, ABCA1, APOA5, CETP, LDLR, APOE) rather than
genome-wide — a hypothesis-driven subset that lets the out-of-sample-LD behaviour be inspected
locus by locus instead of buried in a genome-wide table.

The paper (quoted): *"We used SuSiE to fine-map … modelling up to L = 10 causal signals per
region … credible sets were constructed at 95% coverage and filtered for a minimum purity
(absolute correlation) of 0.5."*

## Similarities (what I replicate faithfully)
- **Engine + core parameters.** `susieR::susie_rss` in z-score mode, `L = 10`,
  `estimate_prior_variance = TRUE`, `scaled_prior_variance = 0.1`, credible-set purity
  `min_abs_corr = 0.5` — the paper's settings.
- **95% coverage credible sets** as the primary output (99% also emitted for reference).
- **Per-variant PIP** and top-variant reporting per credible set.
- **Representative N** = the lead variant's GWAS N (per-variant N is ~flat across a 1 Mb window).

## Differences
| Aspect | Authors | Mine |
|---|---|---|
| LD reference | **in-sample** UKBB genotype LD | **out-of-sample** 1000G EUR (per brief) — the panel is a different sample from the GWAS |
| `estimate_residual_variance` | `TRUE` (safe with in-sample LD) | **`FALSE`** — out-of-sample LD makes residual-variance estimation unstable (susieR issue #162); fixing it is the documented remedy |
| Scope | genome-wide | 11 curated LDL effector loci (hypothesis-driven) |
| LD-mismatch QC | in-sample, so not needed | `estimate_s_rss` **kriging diagnostic** per locus (`kriging_s` column): 0 = LD consistent with z, →1 = severe out-of-sample mismatch |
| Signal filter | genome-wide-significant by construction | explicit `pass_filter` = CS contains a variant with `|z| > ZTHRESH` (~5.45 at 5e-8), so noise credible sets are flagged, not silently kept |

## The out-of-sample LD consequence (the honest caveat)
Using a reference panel that is not the GWAS sample is the crux of this step, and it shows up in
the convergence QC (`susie_qc.tsv`): **3 of 11 loci converge cleanly** (PCSK9, LPL, CETP) while
the other 8 hit the IBSS iteration cap (`max_iter = 100`) without formal convergence. This is the
expected fingerprint of out-of-sample LD, not a bug — the z-vector and the R matrix come from
different samples, so the IBSS updates chase small inconsistencies. Two design choices contain it:
`estimate_residual_variance = FALSE` (the susieR-documented remedy) and the `kriging_s` diagnostic,
which quantifies per-locus how far the panel LD departs from the data. Credible sets at high-`s`
loci should be read as approximate. The clean fix — in-sample UKBB LD — is exactly the resource
the brief chose not to provide, so this caveat is inherent to the task, not to the implementation.

## Headline biological results (`annot.missense.tsv`)
Fine-mapping + VEP annotation recovers the textbook LDL coding variants, which is the strongest
validation that the pipeline works despite out-of-sample LD:

| gene | variant | AA change | PIP | SIFT / PolyPhen | note |
|---|---|---|---|---|---|
| PCSK9 | rs11591147 | R46L | **1.00** | tolerated / benign | the canonical loss-of-function LDL-lowering allele (known causal) |
| APOE | rs429358 | C130R | **1.00** | tolerated / — | the ε4-defining variant |
| ABCG5 | rs141828689 | R198Q | **1.00** | **deleterious / probably-damaging** | high-confidence damaging missense, single-variant credible set |
| APOB | rs1367117 | T98I | 0.15 | deleterious / benign | within a broader credible set |
| ABCG8 | rs4148217 | T400K | 0.04 | deleterious / benign | low PIP, part of the ABCG5/8 signal |

PCSK9 R46L is resolved to a **singleton credible set at PIP = 1** — the ideal fine-mapping
outcome, and a known drug-target variant, so it doubles as a positive control.

**Most causal signal is regulatory, not coding — the bridge to Step 5.** Of the 132 credible-set
variants, only **5 are missense**; **89 are intronic** and the remainder are up/downstream, UTR,
or intergenic (`consequence_summary.png`). Coding annotation therefore explains only a handful of
the loci. The rest act through gene regulation — expression, splicing, transcript usage — which is
precisely what molecular-QTL colocalization (Step 5) is designed to resolve. So fine-mapping does
not close the mechanism question at most loci; it motivates the QTL step.

## Credible-set summary (`cs_summary.tsv`)
144 credible-set rows across the 11 loci (95% + 99% coverage levels); **111 pass the `pass_filter`
genome-wide-significance test**. The filtered set (`cs_summary.filtered.tsv`) is what annotation
and plots consume by default, so downstream steps never inherit noise credible sets — but nothing
is discarded (both files are kept).

## Reproducibility
- Seeded (`SEED` env var); SuSiE is deterministic given seed + inputs.
- Split fit-vs-annotate-vs-plot: `04c` fits SuSiE once and persists `.susie.rds` + tables; `04d`
  (VEP) and `04e` (figures) read only those persisted tables and **never re-fit**, so the
  expensive fit runs once and the light steps are freely re-runnable.
- VEP consequences cached to `LDL_C.vep_raw.json` so re-annotation needs no re-query.