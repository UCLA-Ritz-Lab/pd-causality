# CLAUDE.md — PD Methylation Mendelian Randomization with Causal Foundation Models

## Project goal

Test whether DNA methylation at specific CpGs causally affects Parkinson's disease (PD) risk and severity, using cis-mQTL SNPs as instruments, and test whether tabular foundation models (CausalFM, TabPFN-family PFNs) add anything over classical MR and IV estimators.

Causal question (Mendelian randomization, methylation as exposure):

```
Z (cis-mQTL SNP) ──► X (CpG methylation) ──► Y (PD outcome)
                            ▲                    ▲
                            └──── U (unmeasured confounders) ────┘
```

- **Z (instrument):** cis-mQTL SNPs from the existing PPMI/PEG mQTL results.
- **X (exposure):** blood methylation at the CpG the SNP regulates.
- **Y (outcomes):** PD status; Hoehn & Yahr; MDS-UPDRS (Parts I–III and total); longitudinal UPDRS progression where available.
- **Target:** causal effect of X on Y (and, in Phase 2/3, its heterogeneity), with uncertainty.

Why MR: PD and its treatment (especially levodopa) change blood methylation, so ordinary methylation–outcome associations are uninterpretable. A germline instrument is not affected by disease or treatment.

## IV assumptions to respect (do not skip)

1. **Relevance.** Z must strongly predict X. Report first-stage F-statistics; flag F < 10. Cis-mQTLs are usually strong.
2. **Independence.** Z unconfounded with Y. Adjust for genetic PCs (ancestry), age, sex, and batch/cell composition in the exposure model. Covariates must not be downstream of PD or levodopa: cell composition is fine in the mQTL model but can bias one-sample MR if disease or treatment shifts it.
3. **Exclusion restriction.** Z affects Y only through X. Main threat is pleiotropy: a cis-mQTL may also tag eQTLs and neighbouring CpGs. Use colocalization between the mQTL and PD GWAS signals (`coloc`; coloc-SuSiE where a locus may hold multiple signals). Treat MAPT H1/H2, HLA, GBA1 and LRRK2 regions with special care. Most CpGs have only 1–2 independent cis instruments after clumping; MR-Egger, weighted median and MR-PRESSO need ≥3. Otherwise rely on Wald ratio + coloc, or methods for correlated instruments.
4. **Selection / collider bias.** H&Y and UPDRS exist only in cases, so severity analyses condition on PD status and can open paths between Z and confounders. Progression analyses are subject to index-event bias. State this on every severity result.
5. **Timing.** Record whether methylation was measured before or after treatment. The instrument protects the estimate from reverse causation, but exposure measurement still reflects the timepoint.
6. **First stage in cases.** Case methylation is measured post-diagnosis/treatment, so Z→X can differ by case status. In one-sample analyses, estimate the first stage in controls, or test that it doesn't differ between cases and controls.
7. **Binary outcome scale.** For PD status, state the effect scale (risk difference from 2SLS vs log-OR from ratio methods). ORs are non-collapsible, and case-control sampling distorts risk-difference estimates.

## Where the foundation models fit

- **CausalFM** (Ma et al., ICLR 2026; https://github.com/yccm/CausalFM-toolkit) formalizes back-door, front-door and IV settings, including unobserved confounding. An earlier paper version only instantiated a back-door CATE model, so an IV checkpoint may not be released. Confirm the citation (venue/year) against the repo README in Phase 0.
- **Phase 0 must verify** which settings have released, usable checkpoints, max rows/features, input format, and how uncertainty is reported. Do not assume an IV model is available. If none exists, report that and discuss whether training one is in scope before doing it.
- Justification over 2SLS must be earned empirically: small n (hundreds to low thousands of rows), nonlinear/heterogeneous effects (e.g. by age, sex, PD-relevant variants), Bayesian uncertainty. Report honestly if the PFN does not beat or meaningfully complement classical estimators.

## Workflow

0. **Toolkit recon (done 2026-11-11, see `docs/toolkit_capabilities.md`).** CausalFM cloned to `tools/CausalFM-toolkit` (not pip-installable; add to `sys.path`). torch-CPU + `tabpfn==2.0.9` (pinned, required by CausalFM's fork) in `pd-qtl`. Key findings: only the IV *binary-Z, binary-A* checkpoint works; the continuous-Z checkpoint is broken; there's no continuous-treatment model. `IVModel.estimate_cate` is buggy, so call `.model.estimate_cate_iv(X, A, Y, Z, X_test)`. No GPU, so TabPFN is slow (~1 min per fit).
1. **Nominate (two-sample MR, crossed design).** mQTLs and GWAS must always come from different samples (same-sample pairs are full overlap and bias toward the confounded estimate). Pairings:
   - **Primary:** PPMI and PEG cis-mQTLs (controls preferred) × Nalls 2019 GWAS (GCST009325). Sensitivity: GCST009324 (no proxy cases).
   - **Crossed replication:** PEG mQTLs × PPMI GWAS, and PPMI mQTLs × PEG GWAS.
   - **In-house GWAS** (PPMI and PEG): PLINK2 logistic on imputed dosages, covariates sex, age, ancestry PCs. PPMI genotyped: 511 PD vs 196 HC. Exclude prodromal and SWEDD, and run with and without the 94 PD cases enrolled through LRRK2/GBA/SNCA (ascertained on genotype). These are small; report power and expect wide intervals.
   - Use only full (unfiltered) GWAS summary statistics. A p-filtered file makes missing SNPs non-random and biases MR away from the null. Report how many instruments are missing from the GWAS.
   - Nalls includes PPMI (IPDGC) and possibly PEG samples; note the small overlap.
   - Proxy cases (UK Biobank family history) dilute effects toward the null in GCST009325; compare with GCST009324.
   Filter by first-stage strength and colocalization; control FDR across CpGs. Output: a short list of candidate CpGs with MR estimates, sensitivity analyses, and coloc posteriors.
2. **Benchmark (core contribution).** Semi-synthetic simulation using real genotypes, real methylation and real covariates; simulate Y with known causal effects of X, including nonlinearity, effect heterogeneity, unmeasured X–Y confounding, weak instruments, and pleiotropic Z→Y paths. Confounding can't be added to real X after the fact: either simulate X = f(Z, U, covariates) from real Z and covariates, or define U as a held-out component of real X. Document which. Compare:
   - Classical: 2SLS / Wald ratio, IVW and pleiotropy-robust MR.
   - ML: DML-IV / DeepIV-style estimators.
   - PFN (decided 2026-11-11):
     - **Primary:** TabPFN as the nuisance/stage learner inside a continuous-exposure IV estimator (DML-IV and/or control function), with dosage Z and continuous methylation X.
     - **Comparator:** CausalFM-IV binary checkpoint on dichotomised Z (carrier vs non-carrier) and X (high vs low; state the cut-point). The estimand differs (effect of dichotomised exposure), so report the loss from dichotomising against the same methods run on dichotomised data and on the original scales. Don't treat it as a like-for-like estimate.
   Metrics: bias, RMSE, 95% interval coverage, runtime, behaviour under weak instruments and pleiotropy.
3. **Apply (one-sample, individual-level).** Run the best-performing methods on the Phase 1 shortlist for PD status and progression, in PPMI and PEG separately, then meta-analyse. Selection-bias caveats stated up front. The shortlist was selected using these cohorts' mQTLs, so guard against winner's curse: discover in one cohort and apply in the other (or split-sample). Use baseline / pre-treatment PPMI methylation as the exposure where possible, and state the visit used.

## Data

- **PPMI:** EPIC **v1** array (not v2 — matters for SeSAMe config), repeat-visit methylation, genotypes, clinical outcomes.
- **PEG:** methylation, genotypes, PD status and clinical measures.
- **PD GWAS** (Nalls 2019, full summary statistics from the GWAS Catalog FTP, GRCh37, European, excludes 23andMe), in `data/external/GWAS_PD/`:
  - `GCST009325/GCST009325.tsv`: 15,056 cases + 18,618 proxy cases vs 449,056 controls. Primary.
  - `GCST009324/GCST009324.tsv.gz`: 15,056 cases vs 12,637 controls (no proxies). Sensitivity.
  - GCST009325 columns: `chromosome, base_pair_location, effect_allele, other_allele, beta, standard_error, effect_allele_frequency, p_value, N_cases, N_controls` (beta = log-OR), 17.5M variants.
  - GCST009324 columns: same first 8, then METAL fields (`FreqSE, HetISq, HetChiSq, MaxFreq, HetDf, HetPVal, Direction, MinFreq`); no N columns; **alleles are lowercase** — uppercase before harmonising. 17.4M variants.
  - Neither has an rsID column: match on chr:pos + alleles. Harmonised, tabix-indexed versions (`harmonised/*.h.tsv.gz`) are on the FTP if needed for regional coloc lookups.
  - Plus in-house PPMI/PEG GWAS (see Phase 1).
- **Raw inputs** are symlinked into `data/raw/` from `../PD-QTLs/Downloads/` (PPMI_GCP idats, PPMI_covariates, PPMI_WGS, PEG_450K, PEG_WGS, PEG_covariates, reference).
- **Derived PD-QTLs outputs** are symlinked into `data/pdqtl/` from `../PD-QTLs/results/`. Key files:
  - Imputed genotypes (HRC r1.1, hg19): `PPMI_WGS/post_imputed/ppmi_imputed.{pgen,pvar,psam}` (n=894), `PEG_WGS/post_imputed/peg_imputed.*` (n=1840).
  - PPMI methylation: `PPMI_GCP/local_results/final/beta_matrix.csv.gz` (4.3 GB) + `cell_type_proportions.csv`.
  - cis-mQTLs: `PPMI_MQTL/mqtl/ppmi_{cases,controls}/cis_mqtls.txt.gz`; `PEG_MQTL/mqtl/{peg1_cases,peg1_controls,peg2_cases}/`.
  - Meta-analysis: `META_MQTL/meta/{cases,controls}_meta_results.csv.gz`, `*_sig_fdr05.csv`. Longitudinal: `PPMI_LONG/`.
  - Ancestry PCs: `ANCESTRY_REF/covariates/ancestry_pcs_{PPMI,PEG}.csv`.
- Everything under `data/` is read-only and git-ignored. Write outputs to `results/`, never into the linked trees.
- PPMI IDs: psam IIDs are `PPMISI<PATNO>`; strip the prefix to join to PPMI clinical CSVs (`PATNO`). Case status: `Participant_Status_*.csv` `COHORT_DEFINITION`.
- **PEG methylation is 450K** (`PEG_450K`), PPMI is EPIC v1: restrict to shared probes when combining cohorts.
- Raw and individual-level data never get committed. Paths go in `config/config.yaml`, which is git-ignored; commit `config/config.example.yaml` only.
- No individual-level data or identifiers in logs, notebook outputs, or plots committed to the repo.
- Never send individual-level PPMI/PEG data, or paths to it, to external services, web tools or APIs (PPMI data use agreement).
- Genome builds: EPIC v1 manifest is hg19. Record the build of the mQTL results and the GWAS, and the liftover step if they differ.

## Stack and conventions

- **Pipelines:** Snakemake with per-rule conda environments (`envs/*.yaml`). Default env: `pd-qtl`.
- **R** for statistical genetics and MR: lme4/lmerTest, dplyr, MatrixEQTL, TwoSampleMR / MendelianRandomization, coloc. RMarkdown for documented analyses.
- **Python** for PFN models, simulation, and data wrangling.
- Tools already in use: PLINK2, MatrixEQTL, SeSAMe, missMethyl, bcftools, minimac4.

## Known gotchas

- `readr` silently parses `chr:pos` SNP IDs as times — always specify `col_character()` for ID columns.
- Check zero-variance cell-type covariates on the **genotyped subset**, not the full methylation sample.
- MatrixEQTL pre-filtered output is unsuitable for BH recalculation — use the per-file FDR columns. Note per-file FDR is not a global FDR across chunks; state which FDR is reported.
- Restrict SNPs to the cross-cohort intersection before meta-analysis.
- Harmonise effect alleles and strand between mQTL and GWAS files before any two-sample MR; drop ambiguous palindromic SNPs unless allele frequencies resolve them.
- Handle OOM case by case if it comes up; don't pre-optimise file formats.

## Commands

- `conda activate pd-qtl`
- `snakemake` — run the pipeline (real run; add `-n` only for a dry run).

## Working style

- Debug by reproducing the exact error at the command line on a small input before relaunching a full pipeline.
- Prefer pragmatic, simple solutions; explain any statistical choice in a comment or the RMarkdown text.
- Every result table states: sample, cohort, methylation timepoint, covariates, instrument strength, and which IV assumptions it relies on.

## Suggested layout

```
config/            config.example.yaml (real config.yaml is git-ignored)
data/raw/          symlinks to ../PD-QTLs/Downloads (git-ignored)
data/pdqtl/        symlinks to ../PD-QTLs/results (git-ignored)
data/external/     downloaded external data, e.g. GWAS_PD/ (git-ignored)
envs/              per-rule conda envs
workflow/Snakefile
workflow/rules/    nominate.smk, simulate.smk, estimate.smk, apply.smk
scripts/R/         MR, coloc, mQTL helpers
scripts/py/        PFN wrappers, simulation
analysis/          RMarkdown reports
docs/              toolkit_capabilities.md, methods notes
results/           git-ignored
```
