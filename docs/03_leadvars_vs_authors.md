# Step 3 (lead variants) — my implementation vs. the authors' method

My Step 3: `workflow/Snakefile` rule `lead_variants`, `workflow/scripts/03_leadvars.R`.
Input: the four Step-1 meta tables `results/meta/{trait}.meta.tsv.gz` (~104M variants each).
Output: `results/loci/{trait}.lead.tsv` (independent lead variants) + `{trait}.manhattan.png`.

## What the brief asks vs. what the paper did
The brief lists Step 3 as bare **"Identification of lead variants"** with no LD note (the
"ancestry-appropriate public LD scores" note is attached to Step 4, fine-mapping). So the
brief's Step 3 is **summary-statistics-only**. The paper's fuller method has four sub-steps;
I implement the two that need only summary statistics and document the two that need a
genotype panel as a scoped follow-on.

The paper (quoted): *"For the common and low-frequency variants (MAF > 0.1%, ~15 million
variants), we used the standard genome-wide significance threshold of P < 5×10⁻⁸. However,
our analysis also included up to 80 million rare variants (MAF < 0.1%) for which the standard
P value threshold was too lenient. To be conservative, we treated all rare variant tests as
independent and used the Bonferroni correction to establish a more stringent significance
threshold of P < 0.05/80,000,000 (6.25×10⁻¹⁰). … The variant with the lowest P value was
designated as the lead variant within a 2 Mb locus. In each dataset, neighbouring loci were
merged into one if their lead variants were in LD with an r² of at least 0.05 … we utilized
PLINK to calculate pairwise LD between all lead variants … assigning them into shared
cross-metabolic trait clusters if r² was at least 0.8."*

## Similarities (what I replicate faithfully)
- **Dual frequency-dependent significance threshold.** Common/low-freq (MAF > 0.1%) → P < 5e-8;
  rare (MAF < 0.1%) → P < 6.25e-10 (= 0.05/8e7, the paper's Bonferroni). Config knobs
  `gws_p`, `rare_gws_p`, `maf_rare`.
- **Lowest-P lead within a 2 Mb locus, iterated.** Config `locus_window: 2000000` (±1 Mb).
- **MAF from effect-allele frequency**, taken as min(eaf, 1−eaf) — same frequency the paper
  splits on.

## Differences
| Aspect | Authors | Mine |
|---|---|---|
| Locus independence | distance (2 Mb) **+ r²≥0.05 LD merge** of neighbouring loci | distance (2 Mb) only — no genotype panel provided, so pairwise r² is not computable |
| Cross-trait clustering | PLINK pairwise LD, r²≥0.8 clusters across traits | deferred (needs the same panel) |
| Rare-variant frequency | full MAF coverage from their harmonised panel | fixed the eaf provenance bug first (see below); variants still lacking eaf fall back to the standard line and are flagged `maf_class="unknown"` |
| Thresholding statistic | P value | `neg_log10p` (log-space, so −log₁₀P ≫ 300 loci don't underflow to P=0) |
| Scale handling | cluster | single streaming `zcat \| awk` pass; only significant hits + a thinned Manhattan set enter R (memory-flat on 8 GB) |

## The eaf provenance fix (prerequisite for the dual threshold)
Applying the paper's common-vs-rare split requires a frequency for every variant. In the
original meta output, only EstBB's `eaf` was carried into the join, so the ~76% of variants
seen in UKBB_EUR only had `eaf = NA` and could not be classified. I patched `01_meta.R` to
carry `eaf` from **both** cohorts and combine them (N-weighted mean when both present,
single-cohort value otherwise), mirroring how `beta` is pooled. Step 1 and Step 2 results are
unaffected by this (Step 1 never uses `eaf` in the IVW math; Step 2's per-cohort arms munge
from the raw files, and its meta arm restricts to common HapMap3 SNPs where the affected
rare/single-cohort variants are structurally almost absent) — the incomplete column only bit
at Step 3's threshold, which is why the fix lands here.

## Why distance-only clumping is the honest choice (and its limitation)
The paper's r²≥0.05 merge exists to stop a very strong signal (−log₁₀P ≫ 100) from producing
association that "bleeds" past the fixed 2 Mb window and fragments one true signal into several
apparent loci. Distance clumping cannot detect that two windows tag the same causal variant,
so it can **over-count** loci at the strongest signals. Computing r² needs a matched EUR
genotype panel (raw allele calls for PLINK) — the per-SNP LDSC LD scores in `data/ref/` are a
different object and cannot supply pairwise r². With no panel provided, distance clumping is the
standard, defensible substitute; the r²-merge and cross-trait clustering are an additive
follow-on once a 1000G EUR panel is sourced (also the prerequisite for Step 4 fine-mapping).

## Reproducibility
- Seeded (`params$seed`); the only stochastic element is Manhattan thinning, which is
  deterministic given the seed and does not affect the lead table.
- Sample-tested before deployment: a synthetic fixture exercising common/rare/unknown-eaf
  thresholds, multi-variant loci collapsing to one lead, well-separated loci, and chrX (coded
  23) — all assertions passed (6 leads: 4 common / 1 rare / 1 unknown; the sub-threshold rare
  variant correctly excluded).