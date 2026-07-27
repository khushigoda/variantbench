# VariantBench — metabolite GWAS follow-up pipeline

Config-driven Snakemake + R pipeline reproducing a summary-statistics follow-up of
Tambets et al. for four metabolic traits (LDL_C, total BCAA, glucose, lactate) against two
disease outcomes (CAD, T2D). **European-ancestry throughout** (EstBB + UKBB_EUR exposures,
EUR-matched outcomes, 1000G EUR LD). Exposures are published GWAS summary statistics — no
individual-level data. The four traits form a causal-prior gradient — LDL_C, BCAA, and
glucose as drug-target-validated positive controls, lactate as the weak-prior contrast.

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

**Done (Steps 0–4):**
0. split      — per-chromosome split of each raw `.tsv.gz` (`scripts/00_split_by_chr.sh`);
                streams the meta on an 8 GB box without an OOM. Transient (`temp()`).
1. meta        — fixed-effect IVW meta of EstBB + UKBB_EUR, genome-wide; validate vs published meta
2. ldsc        — SNP-h² + genetic correlation across the four traits (LDSC, Python-2 env)
3. lead_variants — dual-threshold genome-wide-significant loci + Manhattan (LDL_C: 321 loci)
4. finemap     — SuSiE credible sets + VEP missense annotation at 11 LDL effector loci
                 (out-of-sample 1000G EUR LD, kriging QC; PCSK9 R46L & APOE C130R at PIP=1)
5. coloc       — trait–QTL coloc against the eQTL Catalogue, all molecular-QTL quant methods
6. gwmr        — genome-wide MR: metabolite → CAD/T2D (IVW, Egger, median)
7. cismr       — cis (drug-target) MR restricted to effector-gene windows
8. report      — Rmd → self-contained HTML

