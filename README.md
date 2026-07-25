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
workflow/Snakefile        9-step DAG, wildcards expand from config
workflow/scripts/         one R script per step (00–08), snakemake@ object wired
env/*.yml                 conda lockfiles (export the real ones with conda env export)
data/raw/                 downloaded sumstats (git-ignored)
data/ref/                 LD scores, HM3 snplist, 1000G EUR panel (git-ignored)
results/                  all outputs (git-ignored except small tables)
report/report.html        final self-contained report
```

## Run
```bash
conda activate variantbench
snakemake -n                       # dry run — see the whole DAG
snakemake --dag | dot -Tsvg > dag.svg
snakemake --cores 4 --use-conda    # build everything
snakemake --cores 4 results/meta/LDL_C.meta.tsv.gz  # one target (any of LDL_C/BCAA/Glucose/Lactate)
```

## Steps
0. normalize   — raw → common "spine" schema (liftover, −log10p decode, allele harmonize)
1. meta        — fixed-effect IVW meta of EstBB + UKBB_EUR; validate vs published meta
2. ldsc        — h2 + genetic correlation (Python-2 env, Rosetta)
3. lead_variants — genome-wide-significant loci + Manhattan
4. finemap     — SuSiE credible sets + missense/splice annotation (out-of-sample LD!)
5. coloc       — coloc.abf with eQTL Catalogue at a locus
6. gwmr        — genome-wide MR: metabolite → CAD/T2D (IVW, Egger, median)
7. cismr       — cis (drug-target) MR restricted to effector-gene windows
8. report      — Rmd → self-contained HTML

## Key caveat
LD reference is out-of-sample (1000G EUR), not the in-sample UKBB LD the paper used.
Fine-mapping credible sets and clumping can be miscalibrated — stated in the report.
