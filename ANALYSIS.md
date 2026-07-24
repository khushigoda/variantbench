# VariantBench — Analysis Scope & Plan

*Scope locked 2026-07-24. This document defines exactly what is analyzed and why,
before any result exists. Method choices and their rationale are the deliverable;
the numeric results are secondary.*

---

## 1. One-paragraph summary

I take **published GWAS summary statistics** for two circulating metabolites, meta-analyze
two European cohorts, and follow the signal through a standard post-GWAS causal-inference
chain — heritability, locus definition, statistical fine-mapping, colocalization with gene
expression, and Mendelian randomization against two cardiometabolic diseases. No
individual-level genotype data is used. The design is deliberately built around **one
positive-control exposure and one contrast exposure** so the pipeline is shown to
*discriminate*, not just to run.

---

## 2. Exposures — two metabolites

| Metabolite | What it is | Why chosen |
|---|---|---|
| **Total BCAA** | Branched-chain amino acids (valine, leucine, isoleucine), summed | **Positive control.** One of the most reproducible metabolite markers of insulin resistance; elevated circulating BCAAs prospectively predict incident type-2 diabetes. Its genetics are anchored by two well-understood catabolic regulators — **PPM1K** and **BCKDK** — giving a clean effector-gene story that fine-mapping and *cis*-MR should recover. This is the exposure where I expect the whole arc (GWAS → fine-mapped coding gene → coloc → cis-MR → disease) to close. |
| **Lactate** | Glycolytic end-product of anaerobic metabolism | **Contrast.** Different biology (energy metabolism, not amino-acid catabolism), and a weaker/less-established causal link to disease. Its notable locus is **GP6** (platelet collagen receptor GPVI) — biologically distant from a "metabolic" mechanism, which makes it a good test of whether coloc/MR correctly report a *less* compelling causal picture. |

**Why exactly these two (not one, not five):** two exposures on **different chromosomes**
(BCKDK chr16 vs GP6 chr19) and **different pathways** demonstrate the method on a
control-vs-contrast pair while staying inside a 72-hour budget. A single exposure can't
show discrimination; five would spend the time on bookkeeping, not reasoning.

### Locked accessions (GWAS Catalog, GRCh38)

| Metabolite | EstBB | UKBB (EUR) | Published meta (validation) |
|---|---|---|---|
| Total BCAA | `GCST90449544` | `GCST90449793` | `GCST90451287` |
| Lactate    | `GCST90449404` | `GCST90449653` | `GCST90451147` |

The **published meta** column is not an input — it is the target I validate my own
Step-1 meta-analysis against (beta-vs-beta correlation should be ≈1).

---

## 3. Outcomes — two diseases

| Outcome | Source | Type | Role |
|---|---|---|---|
| **T2D** (type-2 diabetes) | Suzuki et al. 2024 (DIAGRAM), EUR | case/control | **Primary** outcome for BCAA — the relationship with the strongest prior evidence. |
| **CAD** (coronary artery disease) | GWAS Catalog `GCST90132314`, EUR-only | case/control | Broader cardiovascular comparator; both metabolites tested against it. |

**Why CAD/T2D:** these are the two cardiometabolic endpoints most strongly tied to BCAA
and lactate biology, both have large well-powered European GWAS with public summary
statistics, and both are case/control so MR returns an interpretable **log-odds-ratio per
SD of metabolite**.

**Deliberate divergence from the reference paper — documented, not accidental:** for CAD I
use the **European-only** study `GCST90132314`, whereas the paper used the multi-ancestry
`GCST90132315`. Reason: my exposure is a European meta-analysis, so an ancestry-matched
outcome keeps LD structure and allele frequencies consistent between exposure and outcome.
Ancestry mismatch is a known source of bias in two-sample MR and LDSC genetic correlation.
This trades some case count for a cleaner causal contrast; the tradeoff is stated in the
report's limitations.

---

## 4. What is fixed vs. open

**Fixed (locked here):** the two metabolites, the two cohorts to meta-analyze, the two
outcomes, the effector genes for *cis*-MR (BCAA → BCKDK, DBT, PPM1K; Lactate → GP6), and
the analysis parameters in `config/config.yaml` (genome-wide p<5e-8, 2 Mb locus window,
MR instrument pruning r²<0.01, ±200 kb cis window, MHC excluded from fine-mapping).

**Open (decided during analysis, once data is in):** the exact fine-mapping loci (driven
by which signals reach significance), the eQTL Catalogue dataset for colocalization (tissue
chosen to match the effector gene), and which MR sensitivity analyses to emphasize.

---

## 5. Report structure — one section per step

The pipeline is nine Snakemake rules (Step 0–8); the report (`workflow/scripts/08_report.Rmd`,
knitting to a single self-contained `report/report.html`) has **one section per analytical
step**, each ending in a plot or table, plus an MCQ appendix.

| § | Report section | Pipeline step | Key output | Biological question it answers |
|---|---|---|---|---|
| 1 | Overview & method choices | — | prose | Why summary-stats; why these exposures/outcomes; the out-of-sample-LD caveat stated up front |
| 2 | Meta-analysis | Step 1 | concordance plot + validation r | Do EstBB and UKBB agree, and does my meta reproduce the published one? |
| 3 | Heritability & genetic correlation | Step 2 (LDSC) | h² ± SE, rg heatmap | How heritable is each metabolite; are BCAA and lactate genetically correlated? |
| 4 | Lead variants | Step 3 | Manhattan + lead table | Where are the genome-wide-significant loci? |
| 5 | Fine-mapping + coding annotation | Step 4 (SuSiE) | credible-set table | Which variant is likely causal, and does it hit a coding effector gene? |
| 6 | Colocalization | Step 5 (coloc) | PP.H4 table + regional plot | Is the metabolite signal driven by the *same* variant as gene expression? |
| 7 | Mendelian randomization | Steps 6–7 | genome-wide MR forest + *cis*-MR table | Does the metabolite *causally* affect CAD/T2D — genome-wide, and via the specific effector gene? |
| 8 | Limitations | — | prose | Out-of-sample LD; single-causal-variant coloc assumption; palindromic SNPs; EUR-only; cis-MR power |
| 9 | MCQ appendix | — | 7 questions | One question per step, each answerable only from the analysis above |

**Why the report reads results rather than computing them:** every heavy computation happens
in its Snakemake rule and writes a small table to `results/`. The report only *reads and
displays* those tables, so it knits in seconds and can't silently recompute something with a
different parameter than the pipeline used. Reproducibility is enforced by the DAG, not by the
notebook.

---

## 6. Headline caveat (carried into every relevant step)

The LD reference for fine-mapping, clumping, and MR pruning is an **out-of-sample** European
panel (1000 Genomes EUR), **not** the in-sample UKBB LD the reference paper used. Out-of-sample
LD can miscalibrate SuSiE credible sets and mis-prune instruments. This is stated in the report
up front and again in each affected step — it is a known limitation I am choosing to accept for
tractability, not an oversight.
