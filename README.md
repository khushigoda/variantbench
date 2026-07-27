# VariantBench — metabolite GWAS follow-up pipeline

Config-driven Snakemake + R pipeline reproducing a summary-statistics follow-up of
Tambets et al. for four metabolic traits (LDL_C, total BCAA, glucose, lactate) against two
disease outcomes (CAD, T2D). **European-ancestry throughout** (EstBB + UKBB_EUR exposures,
EUR-matched outcomes, 1000G EUR LD). Exposures are published GWAS summary statistics — no
individual-level data. The four traits form a causal-prior gradient — LDL_C, BCAA, and
glucose as drug-target-validated positive controls, lactate as the weak-prior contrast.

## Analysis report

**[Full analysis report (HTML)](docs/report.html)** — Plots, tables, fine-mapping credible sets, colocalization results, MR estimates, and MCQ for all steps.

**[Methods summary (ANALYSIS.md)](ANALYSIS.md)** — Scope, exposures, outcomes, and methodological choices (concise overview).

## Layout
```
config/config.yaml        all parameters, accessions, paths  (edit here, not the code)
workflow/Snakefile        DAG,expands from config
workflow/scripts/         R scripts for all 7 steps and related, snakemake@ object wired
data/raw/                 downloaded sumstats (git-ignored)
data/ref/                 LD scores, HM3 snplist, 1000G EUR panel (git-ignored)
results/                  all outputs (git-ignored except small tables)
```

## Status (current submission)

# Steps

**Done:**

0. split       - per-chromosome split of each raw `.tsv.gz`
1. meta        — fixed-effect IVW meta of EstBB + UKBB_EUR, genome-wide;
2. gen corr     — SNP-h² + genetic correlation across the four traits (LDSC)
3. lead_variants — dual-threshold genome-wide-significant loci
4. finemap     — SuSiE credible sets + VEP missense annotation at 11 LDL effector loci
5. coloc       — trait–QTL coloc against the eQTL Catalogue, all molecular-QTL quant methods
6. gwmr        — genome-wide MR: metabolite → CAD/T2D (IVW MR)
7. cismr       — cis (drug-target) MR restricted to effector-gene windows
8. report      — Rmd → self-contained HTML

