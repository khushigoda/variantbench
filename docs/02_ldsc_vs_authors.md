# Step 2 (LDSC) — my implementation vs. the authors' `genetic_correlation` pipeline

Source compared:
https://github.com/ralf-tambets/EstBB-UKBB-metaanalysis/tree/main/code/genetic_correlation
(Nextflow DSL2: `main.nf`, `nextflow.config`, `conf/base.config`, `input/{studies.tsv,metabolites.csv}`)

My Step 2: `workflow/Snakefile` rules `ldsc_munge` / `ldsc_h2` / `ldsc_rg`,
`workflow/scripts/02_ldsc_rg.R`, `scripts/ldsc_prep_hm3.sh`.

## Similarities
- Same tool/version/refs: LDSC 1.0.1, `munge_sumstats.py` -> `ldsc.py --rg`,
  `--merge-alleles w_hm3.snplist`, `--ref-ld-chr`/`--w-ld-chr eur_w_ld_chr/`.
- Staircase rg loop is structurally identical: lead each file over the shrinking tail so
  every C(n,2) pair is computed once. Theirs = main.nf L81-95; mine = 02_ldsc_rg.R L29-36.
- Both parse the "Summary of Genetic Correlation Results" block (they grep `p1 -A N`;
  I readLines + read.table). Same matrix extracted.
- Both feed LDSC the real per-variant N.
- One `--rg` call per lead; LD scores passed as a trailing-slash directory prefix.

## Differences
| Aspect | Authors | Mine |
|---|---|---|
| Orchestrator | Nextflow + SLURM + Singularity | Snakemake + conda, local |
| munge input | per-study regenie GWAS: SNP,A1,A2,FRQ,INFO,N,BETA,SE,P | meta, pre-restricted to SNP,ea,oa,z,n (ldsc_prep_hm3.sh) |
| munge flags | none (default name auto-detect) | explicit --snp/--a1/--a2/--N-col/--signed-sumstats z,0 |
| signed stat | BETA (Z from P-magnitude + BETA-sign) | precomputed z (--signed-sumstats z,0; underflow-safe) |
| QC filters | default --maf-min 0.01 (FRQ) + --info-min 0.9 (INFO) | neither (eaf dropped, no INFO carried); HM3 restriction only |
| heritability | not a deliverable (h2 only inside rg logs) | explicit ldsc_h2 rule + h2.tsv + h2_barplot.png |
| figures | none (text matrix, collectFile) | rg_heatmap.png + h2_barplot.png |
| scope | 249 metabolites x 7 ancestry cohorts | 4 metabolites on the EUR meta (~between_metabolites analog) |

## Notes / open items
- QC gap: authors apply MAF>=0.01 and INFO>=0.9 via munge defaults; I apply neither.
  HM3 restriction keeps ~1.2M common, generally well-imputed variants, so MAF is largely
  redundant, but INFO>=0.9 is a genuine stringency the authors have and I don't.
- P-free munge: my slim file omits pval by design (z fed directly via --signed-sumstats z,0,
  underflow-safe). Unconfirmed without running, because the earlier failed run still had a
  pval column present. If munge errors on a missing p-value, emit pval (col 10) in
  ldsc_prep_hm3.sh — munge auto-detects the name.
