# pd-causality

Does DNA methylation at specific CpGs causally affect Parkinson's disease (PD) risk
and progression? This project uses **Mendelian randomization (MR)**, with germline
cis-mQTL SNPs as instruments for blood methylation, and asks whether tabular
foundation models (TabPFN and CausalFM) add anything over classical IV estimators.

```
Z (cis-mQTL SNP) ──► X (CpG methylation) ──► Y (PD)
                            ▲                   ▲
                            └── U (unmeasured) ─┘
```

PD and its treatment (especially levodopa) change blood methylation, so ordinary
methylation–PD associations can't be interpreted causally. Germline genotype is
fixed at conception and is not changed by disease or treatment.

> **Data access.** The cohort data (PPMI, PEG) are individual-level and under data
> use agreements. None of it is in this repository: no genotypes, methylation,
> phenotypes, identifiers or per-sample results. Only code, configuration templates
> and documentation are committed. Paths are set in a local, git-ignored
> `config/config.yaml`.

## Data

| Source | Role | Notes |
|---|---|---|
| PPMI (EPIC v1 methylation, WGS → HRC r1.1 imputation) | cis-mQTL exposure; later one-sample analyses | baseline methylation: 69 controls, 213 cases (as analysed) |
| PEG (450K methylation, array → HRC r1.1 imputation) | cis-mQTL exposure; later one-sample analyses | PEG1: 231 controls, 312 cases (as analysed); PEG2: 209 cases |
| Nalls et al. 2019 PD GWAS, GWAS Catalog GCST009325 | primary outcome | 15,056 cases + 18,618 UKB proxy cases vs 449,056 controls; GRCh37 |
| Nalls et al. 2019, GCST009324 | sensitivity outcome | no proxy cases (15,056 vs 12,637) |
| 1000 Genomes phase 3 EUR (n = 503) | LD reference for clumping | GRCh37 |

All coordinates are GRCh37/hg19. The upstream mQTL, meta-analysis and longitudinal
mQTL pipelines live in the companion PD-QTLs project.

## Status

### Phase 0: toolkit recon (done)

See [`docs/toolkit_capabilities.md`](docs/toolkit_capabilities.md). CausalFM's only
working IV checkpoint is binary-instrument / binary-treatment, and no
continuous-treatment model exists. The plan is therefore:
- **Primary PFN approach:** TabPFN as the nuisance learner inside a continuous-exposure
  IV estimator (DML-IV / control function).
- **Comparator:** CausalFM-IV on dichotomised Z and X. Its estimand differs, so it is
  not compared like-for-like.

### Upstream QTL fixes (PD-QTLs, done 2026-11)

Before Phase 1 we found and fixed two bugs in the PEG mQTL pipeline:
1. **Sample misalignment.** Genotype columns were in a different order from the
   methylation and covariate columns, so the regression paired different people.
2. **Allele coding.** PEG counted the REF allele, PPMI the ALT allele, so PEG and
   PPMI effects had opposite signs (correlation −0.97).

After the fixes and re-runs, lead mQTL signs agree 100% across cohorts. The PEG/PPMI
meta-analysis, the PPMI longitudinal mQTL analysis and the mQTL×metQTL overlap were
then re-run.

### Phase 1: nominate CpGs by two-sample MR (in progress)

Pipeline in `scripts/R/` (run from the repo root with `conda activate pd-qtl`):

| Step | Script | What it does |
|---|---|---|
| 1 | `01_exposure_candidates.R` | cis-mQTL pairs at p < 5×10⁻⁸ per exposure dataset; drops SNPs within 50 bp of the CpG (probe artefacts); attaches REF/ALT and allele frequency; requires identical alleles across meta-analysed cohorts |
| 2 | `02_harmonise_outcome.R` | matches each SNP to the GWAS on chr:pos + alleles (swaps, strand, palindromic SNPs resolved by frequency or dropped, allele-frequency discordance filter) |
| 3 | `03_clump.R` | greedy LD clumping **separately for each CpG** (r² < 0.001, 1000G EUR) |
| 4 | `04_mr.R` | Wald ratio (1 instrument), IVW with multiplicative random effects (≥ 2), weighted median and MR-Egger (≥ 3), Steiger direction check, Benjamini–Hochberg FDR across CpGs |
| 5 | `05_regional_mqtl.R` | refits **untruncated** regional cis-mQTL statistics (±1 Mb) for the MR hits, for colocalization. Same model as PD-QTLs; reproduces its output to machine precision |
| 6 | `06_coloc.R` | coloc.abf (single causal variant) per hit CpG; p12 = 5×10⁻⁶ (1×10⁻⁵ as sensitivity) |
| 7 | `07_coloc_susie.R` | coloc-SuSiE (several causal variants): PEG1-controls mQTL with in-sample LD, GWAS with 1000G EUR LD |
| 8 | `08_inhouse_gwas.R` | in-house PD GWAS (PLINK2 logistic, European ancestry; sex, age, 5 PCs): PPMI 487 vs 182 (395 vs 182 without genetic-cohort cases), PEG 428 vs 217 |
| 9 | `09_results_table.R` | consolidated table of primary hits: MR, sensitivity, crossed replication, coloc, evidence tier, and the required metadata (sample, timepoint, covariates, instrument strength, IV assumptions) |

Design choices:
- **Exposure.** The primary exposure is the PEG1 + PPMI **controls** mQTL meta-analysis.
  Case methylation is measured after diagnosis and treatment, so the first stage may
  differ in cases (CLAUDE.md, assumption 6). Sensitivity exposures: cases meta,
  PEG1 controls alone, PPMI controls alone.
- **Effect scale.** log-OR of PD per 1.0 change in methylation beta-value. OR per 0.1
  (10 percentage points) is also reported.
- **Instrument availability.** Harmonising before clumping means clumping picks among
  SNPs present in the GWAS. That is an availability criterion, not an outcome-p-value
  one. About 2% of exposure SNPs are missing from or dropped by the GWAS
  harmonisation.

First results (primary: controls meta × GCST009325):
- 12,778 CpGs tested. After clumping, 99.6% have a single independent cis instrument,
  so estimates are almost all Wald ratios. Pleiotropy-robust MR methods need ≥ 3
  instruments and are mostly not applicable; colocalization has to carry that check.
- 44 CpGs pass FDR < 0.05 (27 Bonferroni), grouped into 14 loci.
  - Half of the loci are known Nalls 2019 PD risk loci, including the MAPT H1/H2
    region, which is flagged for special care.
  - All 44 keep their direction and p < 0.05 in the no-proxy GWAS.
  - They also keep their direction with case-derived and single-cohort mQTLs.
- Genomic inflation of the per-CpG MR statistics: λ ≈ 1.2.

Colocalization. A Wald ratio from one SNP can't tell a shared causal variant from
separate signals in LD, so each hit is tiered by coloc:

| Tier | Meaning | CpGs | Loci |
|---|---|---|---|
| A | shared variant: coloc-SuSiE H4 ≥ 0.8, or coloc.abf H4 ≥ 0.8 where SuSiE had no credible set | 5 | 5 |
| B | moderate shared-variant support (H4 0.5–0.8) | 11 | 4 |
| C | distinct signals (H3 favoured), i.e. likely LD confounding | 28 | 5 |

- The largest known PD loci among the hits (including MAPT and SNCA) fall in tier C,
  even with multi-signal SuSiE.
- Several tier A/B loci rest on sub-genome-wide PD signals, where SuSiE finds no
  GWAS credible set. They depend on coloc.abf and need replication.

Crossed replication with in-house GWAS (independent of Nalls). The in-house GWAS have
fewer than 500 cases, so they test direction, not significance:
- **PPMI mQTL × PEG GWAS:** 92% of CpGs with Nalls-based MR p < 1e-5 have the same
  direction (cor(z) = 0.78).
- **PEG mQTL × PPMI GWAS:** 77–83% same direction (cor(z) ≈ 0.7) once the MAPT region
  is excluded. PPMI's GWAS shows essentially no MAPT effect (lead-SNP log-OR 0.006 vs
  0.24 in Nalls), which reverses the concordance when MAPT is included.

These are **nominations, not causal claims**. Remaining Phase 1 work:
- [x] Colocalization (coloc.abf, coloc-SuSiE).
- [x] Crossed replication with in-house GWAS.
- [x] Results table with sample, cohort, methylation timepoint, covariates, instrument
      strength and the IV assumptions relied on (`results/phase1/tables/`).
- [ ] CpGs sharing one instrument: MR alone cannot tell which CpG, if any, is causal.
- [ ] MAPT H1/H2 handling (inversion; extended LD).
- [ ] Coloc-SuSiE for exposures other than PEG1 controls, and for the no-proxy GWAS.
- [ ] Larger PEG GWAS if a phenotype file for all ~1,840 genotyped PEG participants
      becomes available (currently only the 760 methylation-study participants).

### Phase 2: semi-synthetic benchmark (planned)

Real genotypes, methylation and covariates, with simulated outcomes. The simulations
cover known nonlinear and heterogeneous effects, unmeasured confounding, weak
instruments and pleiotropy. Methods compared: classical MR / 2SLS, DML-IV / DeepIV,
TabPFN-in-IV, and CausalFM-IV (dichotomised).

### Phase 3: one-sample application (planned)

Run the best methods on the Phase 1 shortlist in PPMI and PEG separately (status and
progression), then meta-analyse. To avoid winner's curse, discover in one cohort and
apply in the other.

## Layout

```
config/config.example.yaml   paths + thresholds (copy to config.yaml, git-ignored)
scripts/R/                   Phase 1 MR pipeline (phase1_lib.R + numbered steps)
scripts/py/                  PFN / CausalFM wrappers (Phase 0 smoke test)
docs/                        toolkit capabilities, methods notes
data/                        symlinks to raw + PD-QTLs outputs (git-ignored)
results/                     all outputs (git-ignored)
```

See [`CLAUDE.md`](CLAUDE.md) for the full analysis plan, the IV assumptions checked
at each step, and the known data gotchas.
