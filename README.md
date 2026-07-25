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
workflow/Snakefile        DAG (Steps 0-2 active; 3-8 deferred, see below), expands from config
workflow/scripts/         one R script per step (00–08), snakemake@ object wired
env/*.yml                 conda lockfiles (export the real ones with conda env export)
data/raw/                 downloaded sumstats (git-ignored)
data/ref/                 LD scores, HM3 snplist, 1000G EUR panel (git-ignored)
results/                  all outputs (git-ignored except small tables)
report/report.html        final self-contained report
```

## Status (current submission)

**Step 1 — meta-analysis: validated.** The from-scratch fixed-effect IVW meta reproduces the
authors' published `meta_EUR` essentially exactly. On **BCAA** (first trait run end-to-end,
genome-wide, all 103,927,003 union variants):

| metric | value |
|---|---|
| Pearson r (our β vs published β) | **0.999999999999989** (≈1.0000) |
| OLS slope | **0.999999999987551** (≈1.0000) |
| variants compared to published | 103,927,003 |
| in both cohorts / EstBB-only / UKBB-only | 17,976,457 / 8,197,137 / 77,753,409 |

Outputs in `results/meta/BCAA.{meta.tsv.gz, validation.tsv, validation.png, concordance.png}`.
The remaining three traits (Lactate, LDL_C, Glucose) run with the identical rule and are being
produced by the per-trait loop below. **Step 2 (LDSC h²/rg)** scripts are prepared and unit-tested;
they run once the four metas exist.

## Run

Per-step conda envs are activated by Snakemake itself via `--use-conda` (R steps →
`variantbench-r`; LDSC → `ldsc`), so no manual `conda activate` is needed. Run from the repo root.

```bash
snakemake -n                          # dry run — see the whole DAG
snakemake --dag | dot -Tsvg > dag.svg # draw it

# --- Step 1 (meta), one trait at a time (8 GB / ~30 GB-disk safe) ---
# --cores 1 keeps peak memory to a single chromosome; targeting ONE trait at a time keeps
# peak disk to that trait's transient splits (temp() deletes them once its meta succeeds).
for T in BCAA Lactate LDL_C Glucose; do
  snakemake --use-conda --cores 1 results/meta/${T}.meta.tsv.gz || break
done

# --- Step 2 (LDSC h2 + rg) once all four metas exist ---
snakemake --use-conda --cores 1 results/ldsc/rg.tsv results/ldsc/h2.tsv \
          results/ldsc/rg_heatmap.png results/ldsc/h2_barplot.png

# Or the whole thing (Steps 0-2) via the default target — but on a memory/disk-constrained
# box prefer the per-trait loop above so splits don't pile up across traits:
# snakemake --use-conda --cores 1
```

## Steps

**Active now (Steps 0–2 — the current submission scope):**
0. split      — per-chromosome split of each raw `.tsv.gz` (`scripts/00_split_by_chr.sh`);
                streams the meta on an 8 GB box without an OOM. Transient (`temp()`).
1. meta        — fixed-effect IVW meta of EstBB + UKBB_EUR, genome-wide; validate vs published meta
2. ldsc        — SNP-h² + genetic correlation across the four traits (LDSC, Python-2 env)

**Deferred (Steps 3–8 — scaffolded on disk, rules commented out at the foot of the Snakefile;
restore by re-enabling the disabled `config` keys, see the header comment there):**
3. lead_variants — genome-wide-significant loci + Manhattan
4. finemap     — SuSiE credible sets + missense/splice annotation (out-of-sample LD!)
5. coloc       — coloc.abf with eQTL Catalogue at a locus
6. gwmr        — genome-wide MR: metabolite → CAD/T2D (IVW, Egger, median)
7. cismr       — cis (drug-target) MR restricted to effector-gene windows
8. report      — Rmd → self-contained HTML

*Note: `scripts/00_normalize.R` (outcome harmonization + T2D GRCh37→GRCh38 liftover) belongs to
the deferred outcome path, not the active metabolite path — its Step-0 role there is distinct
from the per-chromosome split used for the exposures above.*

## Key caveat
LD reference is out-of-sample (1000G EUR), not the in-sample UKBB LD the paper used.
Fine-mapping credible sets and clumping can be miscalibrated — stated in the report.
