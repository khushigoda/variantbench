# VariantBench — Analysis Scope & Plan

*Scope locked 2026-07-24. This document defines exactly what is analyzed and why,
before any result exists. Method choices and their rationale are the deliverable;
the numeric results are secondary.*

---

## 1. One-paragraph summary

I take **published GWAS summary statistics** for four circulating metabolic traits, meta-analyze
two European cohorts, and follow the signal through a standard post-GWAS causal-inference
chain — heritability, locus definition, statistical fine-mapping, colocalization with gene
expression, and Mendelian randomization against two cardiometabolic diseases. No
individual-level genotype data is used, and **the analysis is European-ancestry throughout**
(exposures, disease outcomes, and LD reference are all EUR) so LD structure and allele
frequencies stay consistent across every two-sample step. The design is deliberately built
around a **gradient of causal-prior strength** — three exposures with well-established,
drug-target-validated links to disease as positive controls, and one contrast exposure — so
the pipeline is shown to *discriminate*, not just to run.

---

## 2. Exposures — four metabolic traits (EUR only)

| Trait | What it is | Role & why chosen |
|---|---|---|
| **LDL_C** | LDL cholesterol (measured, not clinically calculated) | **Positive control — the canonical MR benchmark.** The causal LDL→CAD relationship is proven by RCTs of statins (HMGCR), PCSK9 inhibitors, and ezetimibe (NPC1L1). If the pipeline recovers a positive LDL_C→CAD *cis*-MR at these drug-target genes, that is the strongest possible validation the whole arc works — every effector gene is a real therapeutic target. |
| **Total BCAA** | Branched-chain amino acids (valine, leucine, isoleucine), summed | **Positive control.** One of the most reproducible metabolite markers of insulin resistance; elevated circulating BCAAs prospectively predict incident type-2 diabetes. Its genetics are anchored by the catabolic regulators **PPM1K** and **BCKDK**, giving a clean effector-gene story fine-mapping and *cis*-MR should recover — this mirrors the reference email's GOOD MCQ example. |
| **Glucose** | Fasting circulating glucose | **Positive control.** Glycemia is mechanistically central to T2D; fasting-glucose genetics are anchored by glucokinase (**GCK**), its regulator **GCKR**, islet **G6PC2**, and the ZnT8 transporter **SLC30A8**. Glucose→T2D should recover a positive effect, and **GCK is itself a drug target** (glucokinase activators), keeping the drug-target *cis*-MR interpretation intact. |
| **Lactate** | Glycolytic end-product of anaerobic metabolism | **Contrast.** Different biology (energy metabolism, not lipid / amino-acid / glycemic), a weaker causal prior, and a notable locus at **GP6** (platelet collagen receptor GPVI) that is biologically distant from a metabolic disease mechanism — a test of whether coloc/MR correctly report a *less* compelling causal picture. |

**Why these four (a gradient, not near-duplicates):** the exposures sit on **four different
chromosomes** and **four different pathways** (lipid, amino-acid catabolism, glycemia,
glycolysis), spanning strong-and-RCT-validated (LDL_C) → strong (BCAA, Glucose) → weak
(Lactate) causal priors. Recovering the three controls *and* correctly down-weighting the
contrast is what demonstrates discrimination; a single exposure can't show it, and
near-duplicate traits would only re-test the same locus.

### Locked accessions (GWAS Catalog, GRCh38 — EUR cohorts only)

| Trait | EstBB | UKBB (EUR) | Published meta (validation) |
|---|---|---|---|
| LDL_C      | `GCST90449408` | `GCST90449657` | `GCST90451151` |
| Total BCAA | `GCST90449544` | `GCST90449793` | `GCST90451287` |
| Glucose    | `GCST90449379` | `GCST90449628` | `GCST90451122` |
| Lactate    | `GCST90449404` | `GCST90449653` | `GCST90451147` |

All from Supplementary Table S12. The 6 non-EUR UKBB ancestries (AFR/AMR/CSA/EAS/MID) and the
all-ancestry meta are **deliberately not used** — the analysis is EUR-only. The **published
meta** column is not an input; it is the target I validate my own Step-1 meta-analysis against
(beta-vs-beta correlation should be ≈1). LDL_C is the S12 trait "LDL cholesterol", *not* the
separate "Clinical LDL cholesterol".

---

## 3. Outcomes — two diseases

| Outcome | Source | Type | Role |
|---|---|---|---|
| **CAD** (coronary artery disease) | GWAS Catalog `GCST90132314` (EUR-only, primary) + `GCST90132315` (multi-ancestry, sensitivity), harmonised GRCh38 | case/control | **Primary** outcome for LDL_C — the RCT-proven relationship. Cardiovascular comparator for the other traits. |
| **T2D** (type-2 diabetes) | Suzuki et al. 2024 (DIAGRAM), EUR (primary) + all-ancestry (sensitivity) | case/control | **Primary** outcome for BCAA and Glucose — the diabetes-linked exposures. |

**Why CAD/T2D:** these two cardiometabolic endpoints map cleanly onto the four exposures
(LDL_C→CAD; BCAA/Glucose→T2D; Lactate weakly to both), both have large well-powered
**European** GWAS with public summary statistics, and both are case/control so MR returns an
interpretable **log-odds-ratio per SD of exposure**. Keeping the outcomes European matches the
EUR exposures and avoids ancestry-mismatch bias in two-sample MR.

**Deliberate choice on the CAD accession — documented, not accidental:** the paper's Methods
cite the multi-ancestry study `GCST90132315`, but the authors' *released code* actually reads
the **European-only** `GCST90132314` (their `Aragam_2022_GCST90132314_harmonized.parquet`) —
a paper-vs-code discrepancy I verified against their GitHub repo. I pull **both**: I use the
EUR-only `GCST90132314` as the **primary** outcome (ancestry-matched to my European exposure,
keeping LD structure and allele frequencies consistent — ancestry mismatch is a known source
of bias in two-sample MR and LDSC), and keep the multi-ancestry `GCST90132315` as a
**sensitivity check** (more cases, larger power). The tradeoff is stated in the report's
limitations.

**Genome build / liftover:** I download the GWAS Catalog **harmonised** CAD files, which are
already lifted to GRCh38 (`genome_assembly: GRCh38` in the file metadata), matching my GRCh38
metabolite exposures — so **CAD needs no liftover**. T2D is different: DIAGRAM serves the
Suzuki 2024 sumstats in **GRCh37/hg19 only**, so T2D **still requires a GRCh37→GRCh38 liftover**
in Step 0 before MR harmonisation. This asymmetry is intentional and handled by the per-outcome
`build:` field in the config (Step 0 lifts only when `build == GRCh37`).

---

## 4. What is fixed vs. open

**Fixed (locked here):** the four metabolic traits, the two EUR cohorts to meta-analyze, the two
outcomes, the effector genes for *cis*-MR (LDL_C → PCSK9, HMGCR, NPC1L1, LDLR; BCAA → BCKDK,
DBT, PPM1K; Glucose → GCK, GCKR, G6PC2, SLC30A8; Lactate → GP6), and the analysis parameters
in `config/config.yaml` (genome-wide p<5e-8, 2 Mb locus window, MR instrument pruning r²<0.01,
±200 kb cis window, MHC excluded from fine-mapping).

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
| 3 | Heritability & genetic correlation | Step 2 (LDSC) | h² ± SE, rg heatmap | How heritable is each trait; how are the four traits genetically correlated (e.g. Glucose–BCAA insulin-resistance axis vs Lactate)? |
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
