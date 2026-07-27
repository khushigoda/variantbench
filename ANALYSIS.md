# VariantBench — Methods

*Analysis pipeline executing steps 0–8: meta-analysis, heritability, fine-mapping, colocalization, and MR. Scope locked 2026-07-24.*

---

## 1. Overview

Four circulating metabolic traits (LDL_C, BCAA, Glucose, Lactate) — meta-analyzed across EstBB and UKBB-EUR cohorts — are interrogated through LDSC, SuSiE fine-mapping, eQTL colocalization, and two-sample MR against CAD and T2D. Analysis is EUR-only throughout (ancestry-matched exposures, outcomes, and 1000G-EUR LD reference). Design includes three drug-target positive controls (LDL_C, BCAA, Glucose) and one contrast (Lactate) to demonstrate discriminative power.

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

## 5. Pipeline structure

The pipeline is nine Snakemake rules (Step 0–8); the report (`workflow/scripts/08_report.Rmd`,
knitting to a single self-contained `report/report.html`) has **one section per analytical
step**, each ending in a plot or table, plus an MCQ appendix.

| § | Report section | Pipeline step | Key output | Status |
|---|---|---|---|---|
| 1 | Overview & method choices | — | prose | ✓ Completed |
| 2 | Meta-analysis | Step 1 | concordance plot + validation r | ✓ Completed (EstBB + UKBB-EUR → meta.tsv) |
| 3 | Heritability & genetic correlation | Step 2 (LDSC) | h² ± SE, rg heatmap | ✓ Completed (h²_observed, rg matrix) |
| 4 | Lead variants | Step 3 | Manhattan + lead table | ✓ Completed (321 leads, distance-clumped) |
| 5 | Fine-mapping + coding annotation | Step 4 (SuSiE) | credible-set table | ✓ Completed (L=10, coverage=0.95, 1000G-EUR LD) |
| 6 | Colocalization | Step 5 (coloc) | PP.H4 table + regional plot | ✓ Completed (eQTL Catalogue, tissue-matched) |
| 7 | Mendelian randomization | Steps 6–7 | genome-wide MR forest + *cis*-MR table | ✓ Completed (IVW + cis-MR drug-target validation) |
| 8 | Limitations | — | prose | ✓ Completed (out-of-sample LD, power, ancestry) |
| 9 | MCQ appendix | — | 7 questions | ✓ Completed (one per analytical step) |

**Why the report reads results rather than computing them:** every heavy computation happens
in its Snakemake rule and writes a small table to `results/`. The report only *reads and
displays* those tables, so it knits in seconds and can't silently recompute something with a
different parameter than the pipeline used. Reproducibility is enforced by the DAG, not by the
notebook.

---

## 6. Key methodological choices

| Aspect | Choice | Rationale |
|---|---|---|
| **LD reference** | 1000G Phase 3 EUR (526 founders, hg38), out-of-sample | Consistent ancestry across exposures, outcomes, and LD. Out-of-sample is a known limitation (miscalibration risk), accepted for tractability. estimate_s_rss diagnostic flags severity. |
| **Ancestry** | EUR-only throughout | Exposures (EstBB+UKBB-EUR), outcomes (CAD EUR, T2D EUR), LD (1000G-EUR). Avoids ancestry-mismatch bias in MR and LDSC. |
| **Fine-mapping** | SuSiE on 1000G-EUR LD: L=10, coverage=0.95, scaled_prior_variance=0.1, estimate_residual_variance=FALSE | Matches authors' (Tambets et al.) recommendations; out-of-sample LD noted as limitation. |
| **MR instruments** | Clumped to r²<0.01 against 1000G-EUR; palindromes excluded; harmonized on chr:pos+allele | Standard practice; out-of-sample LD may prune differently than in-sample. |
| **Positive controls** | LDL_C, BCAA, Glucose (drug-target validated) vs. Lactate (contrast) | Tests pipeline discrimination: recovers strong signal at known targets; down-weights weaker biology. |
| **Colocalization** | eQTL Catalogue, tissue-matched | Shared-variant assumption within single eQTL per locus; tissue choice matches effector gene biology. |

---

## 7. Rendered report

Complete analysis (plots, tables, MCQ): **[report.html](docs/report.html)**
