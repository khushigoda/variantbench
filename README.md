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

## Viewing the report (no rerun needed)

`report/report.html` is **fully self-contained** — every figure is base64-inlined and all
JS/CSS is embedded, so it needs no internet and no pipeline run to view. GitHub's file
browser shows HTML *source* rather than rendering it, so use one of:

- **One click (in-browser):**
  [**View the report via htmlpreview**](https://htmlpreview.github.io/?https://github.com/khushigoda/variantbench/blob/main/report/report.html)
  *(the file is ~6.7 MB, so give it a few seconds to load).*
- **Always works (offline):** clone or download the repo and open `report/report.html`
  in any browser — e.g. `open report/report.html` (macOS) / `xdg-open report/report.html` (Linux).

Regenerate from source at any time with `snakemake --use-conda --cores 1 --force report`.

## Status (current submission)

**Step 1 — meta-analysis: validated for all four traits.** The from-scratch fixed-effect IVW
meta reproduces the authors' published `meta_EUR` at machine precision, genome-wide (all
103,927,003 union variants per trait — every variant kept, not HapMap3-restricted):

| trait | Pearson r (our β vs published β) | OLS slope |
|---|---|---|
| BCAA    | 0.999999999999989 | 0.999999999987551 |
| Lactate | 0.999999999999991 | 1.00000000001978 |
| LDL_C   | 0.999999999999986 | 1.00000000002351 |
| Glucose | 1.000000000000000 | 1.00000000000254 |

Every trait draws on the same 103,927,003-variant union (17,976,457 in both cohorts;
8,197,137 EstBB-only; 77,753,409 UKBB-only) — EstBB and UKBB_EUR release the same imputed
variant set per trait, so only the effect sizes differ. Outputs per trait in
`results/meta/{trait}.{meta.tsv.gz, validation.tsv, validation.png, concordance.png}`.

**Step 2 (LDSC h²/rg)** replicates the paper's two genetic-correlation analyses over **three
datasets × four traits**: the two raw per-cohort GWAS (`EstBB`, `UKBB_EUR`) and our own Step-1
meta (`meta_EUR`, which carries native rsID). It produces:

- **Analysis A — rg between biobanks, per trait:** for each trait, correlate the 3 datasets
  → C(3,2)=3 × 4 = **12 rg** (`rg_between_biobanks.tsv`).
- **Analysis B — rg between traits, per dataset:** within each dataset, correlate the 4 traits
  → C(4,2)=6 × 3 = **18 rg** (`rg_between_traits.tsv`).
- **SNP-h²** for all 12 (trait, dataset) pairs (`h2.tsv`), plus three figures.

**QC is munge's own**, matching the authors (who pass no filter flags). `--merge-alleles
w_hm3.snplist` does the HapMap3 restriction + allele harmonisation; munge's defaults
`--maf-min 0.01` and `--info-min 0.9` fire automatically from the `FRQ`/`INFO` columns. INFO≥0.9
therefore applies to the two raw cohorts (which carry `info`); `meta_EUR` has no INFO column, so
it gets HM3+MAF only — exactly how INFO was handled for the published meta (a per-cohort upstream
QC, never a single meta-stage filter). A thin adapter (`scripts/ldsc_munge_input.sh`, **no
MAF/INFO filtering**) is the only pre-step: it converts `neg_log10(P)` → real P (1e-300 floor),
writes `n` as a float (the LDSC env's Python-2.7 pandas rejects an integer N column), renames
columns to canonical `SNP/A1/A2/BETA/P/N/FRQ/INFO`, **and pre-restricts to the ~1.2M HapMap3
rsIDs** via a streaming `awk` join against `w_hm3.snplist`. That last step is a performance
fix, not a science change: munge's Py2.7 pandas would otherwise parse all 26–104M genome-wide
rows before `--merge-alleles` discards the ~99% non-HapMap3 ones.

Two Step-2 parameters keep the run tractable on an 8 GB laptop. Neither changes the result —
`--merge-alleles` still performs all allele harmonisation inside munge and the FRQ/INFO QC
defaults still fire:
1. **HapMap3 pre-restriction (adapter):** bounds the file munge parses to the ~1.2M HM3 rsIDs,
   so munge never materialises the 26–104M genome-wide rows only to discard 99% of them.
2. **`--chunksize 200000` on the munge call:** munge's default chunksize (`5e6`) exceeds the
   ~1.2M-row HM3 file, which forces the `--merge-alleles` membership test to run once over the
   whole file. Under LDSC's pinned pandas 0.20.3 that single large-vs-large string comparison is
   super-linear; processing the file in 200k-row chunks makes the same comparison run in small
   batches (`7 × isin(200k) ≪ 1 × isin(1.2M)`), reducing each munge job to ~20 s. A chunksize
   above 1.2M has no effect (still one chunk) — 200k is the load-bearing value.

Signed statistic = `BETA` (matches the authors).

**Step 3 — lead variants: done for LDL_C.** Dual frequency-dependent genome-wide-significance
threshold (common/low-freq MAF>0.1% → P<5e-8; rare → P<6.25e-10, the paper's Bonferroni),
lowest-P lead within a 2 Mb locus. **321 independent lead loci** on the LDL meta, plus a
Manhattan figure. Distance-only clumping (no genotype panel provided for the paper's r²≥0.05
merge). Outputs: `results/loci/LDL_C.{lead.tsv, manhattan.png}`. Method vs. authors:
`docs/03_leadvars_vs_authors.md`.

**Step 4 — fine-mapping: done for the 11 LDL effector loci.** `susie_rss` z-mode, L=10,
`estimate_residual_variance=FALSE` (out-of-sample 1000G EUR LD, per the brief),
purity `min_abs_corr=0.5`, 95%/99% coverage, with an `estimate_s_rss` (kriging) LD-mismatch
diagnostic per locus. **144 credible sets (111 pass the genome-wide-sig filter)** across
PCSK9/APOB/ABCG5-8/HMGCR/NPC1L1/LPL/ABCA1/APOA5/CETP/LDLR/APOE. VEP annotation recovers the
textbook coding variants — **PCSK9 R46L and APOE C130R each at PIP=1**, ABCG5 R198Q (damaging)
at PIP=1. Outputs: `results/finemap/LDL_C.{finemap.tsv, cs_summary[.filtered].tsv, susie_qc.tsv,
annot[.missense].tsv}` + 11 regional plots and two summary figures. Method vs. authors:
`docs/04_finemap_vs_authors.md`.

**Step 5 — colocalization with molecular QTLs: in design.** Trait–QTL coloc of the LDL
fine-mapping signals against the eQTL Catalogue, sweeping all molecular-QTL quantification
methods (ge/exon/tx/txrev/leafcutter/microarray/aptamer), tissue chosen data-driven per gene.
Retrieval + engine architecture under review (see `docs/05_coloc_design.md`).

## Run

Three conda envs, nested: you activate only the **driver** env (`variantbench`, which holds
`snakemake`); Snakemake then dispatches each rule into its own env. R steps run in
`variantbench-r` (activated by `--use-conda`); the LDSC binary runs in `ldsc` (Python-2.7),
reached by explicit `conda run -n ldsc` inside the Step-2 rules. So the only `conda activate`
you type is the driver; `--use-conda` handles the rest. Run from the repo root.

```bash
conda activate variantbench           # the driver env — the one with snakemake
```

```bash
snakemake -n                          # dry run — see the whole DAG
snakemake --dag | dot -Tsvg > dag.svg # draw it

# --- Step 1 (meta), one trait at a time (8 GB / ~30 GB-disk safe) ---
# --cores 1 keeps peak memory to a single chromosome; targeting ONE trait at a time keeps
# peak disk to that trait's transient splits (temp() deletes them once its meta succeeds).
for T in BCAA Lactate LDL_C Glucose; do
  snakemake --use-conda --cores 1 results/meta/${T}.meta.tsv.gz || break
done

# --- Step 2 (LDSC h2 + both rg analyses) once all four metas exist ---
# 12 munge jobs (3 datasets x 4 traits) + --h2 + the two rg staircases. The adapter pre-restricts
# each genome-wide input (26M-104M variants) to the ~1.2M HapMap3 rsIDs before munge, and the rule
# passes --chunksize 200000 so munge's pinned-pandas allele-merge runs in small batches. Together:
# ~20s/munge job, a few hundred MB RAM, and --merge-alleles + FRQ/INFO QC unchanged inside munge.
# --rerun-triggers mtime: trust the existing frozen metas (skip re-running Step 1); it MUST go
# LAST because it consumes multiple values and would otherwise swallow the target paths.
snakemake --use-conda --cores 6 \
          results/ldsc/h2.tsv \
          results/ldsc/rg_between_biobanks.tsv results/ldsc/rg_between_traits.tsv \
          results/ldsc/h2_barplot.png \
          results/ldsc/rg_between_biobanks.png results/ldsc/rg_between_traits.png \
          --rerun-triggers mtime

# Or the whole thing (Steps 0-2) via the default target — but on a memory/disk-constrained
# box prefer the per-trait loop above so splits don't pile up across traits:
# snakemake --use-conda --cores 1
```

## Steps

**Done (Steps 0–4):**
0. split      — per-chromosome split of each raw `.tsv.gz` (`scripts/00_split_by_chr.sh`);
                streams the meta on an 8 GB box without an OOM. Transient (`temp()`).
1. meta        — fixed-effect IVW meta of EstBB + UKBB_EUR, genome-wide; validate vs published meta
2. ldsc        — SNP-h² + genetic correlation across the four traits (LDSC, Python-2 env)
3. lead_variants — dual-threshold genome-wide-significant loci + Manhattan (LDL_C: 321 loci)
4. finemap     — SuSiE credible sets + VEP missense annotation at 11 LDL effector loci
                 (out-of-sample 1000G EUR LD, kriging QC; PCSK9 R46L & APOE C130R at PIP=1)

**In design (Step 5) / scaffolded (Steps 6–8 — rules commented at the foot of the Snakefile):**
5. coloc       — trait–QTL coloc against the eQTL Catalogue, all molecular-QTL quant methods
6. gwmr        — genome-wide MR: metabolite → CAD/T2D (IVW, Egger, median)
7. cismr       — cis (drug-target) MR restricted to effector-gene windows
8. report      — Rmd → self-contained HTML

*Note: `scripts/00_normalize.R` (outcome harmonization + T2D GRCh37→GRCh38 liftover) belongs to
the deferred outcome path, not the active metabolite path — its Step-0 role there is distinct
from the per-chromosome split used for the exposures above.*

## Key caveat
LD reference is out-of-sample (1000G EUR), not the in-sample UKBB LD the paper used.
Fine-mapping credible sets and clumping can be miscalibrated — stated in the report.
