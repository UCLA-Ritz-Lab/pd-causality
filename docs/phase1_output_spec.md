# Phase 1 output specification

Column-by-column description of every file Phase 1 writes under `results/phase1/`
(git-ignored). Files are tab-separated; `.tsv.gz` files are gzip-compressed.

## Conventions used in all files

| Item | Meaning |
|---|---|
| Coordinates | GRCh37/hg19 throughout. |
| `snp_id` | `chr:pos` (no alleles). Always match on alleles too: see `ref` / `alt`. |
| `cpg_id` | Illumina probe ID (450K / EPIC v1). |
| Effect allele | **ALT** (from the HRC-imputed `.pvar`) for every mQTL and harmonised GWAS effect. `eaf` = ALT allele frequency. |
| mQTL `beta` | Change in methylation **beta-value** (0–1 scale) per ALT allele. |
| GWAS `beta_y` | log-OR of PD per ALT allele. |
| MR `b` | log-OR of PD per 1.0 change in methylation beta-value. `or_per_0.1` = exp(0.1·b) is the OR per 10 percentage points. |
| `{exposure}__{outcome}` | File-name pattern for an analysis ("combo"), e.g. `meta_controls__nalls2019`. |
| `flag_region` | Long-range-LD / pleiotropy-prone region the SNP or CpG lies in: `MAPT_17q21.31` (17:43.4–44.9 Mb), `HLA`, `GBA1` (1:154.2–156.2 Mb), `LRRK2` (12:39.6–41.8 Mb); empty/NA otherwise. |
| `fdr` | Benjamini–Hochberg across CpGs within one exposure × outcome analysis. |

**Exposures (mQTL datasets)**

| Name | Content |
|---|---|
| `meta_controls` | **Primary.** Common-effect IVW meta of PEG1 controls (n = 231, 450K) + PPMI controls (n = 69, EPIC v1), from PD-QTLs `META_MQTL`. |
| `meta_cases` | Same for cases (PEG1 312 + PPMI 213). Sensitivity: case methylation is post-diagnosis. |
| `peg1_controls`, `ppmi_controls` | Single-cohort controls mQTLs (PD-QTLs MatrixEQTL output). |

**Outcomes (PD GWAS)**

| Name | Content |
|---|---|
| `nalls2019` | **Primary.** GCST009325: 15,056 cases + 18,618 UKB proxy cases vs 449,056 controls. |
| `nalls2019_noproxy` | GCST009324: 15,056 cases vs 12,637 controls (no proxies). |
| `ppmi_gwas`, `ppmi_nogen_gwas` | In-house PPMI GWAS (step 8), with / without LRRK2/GBA/SNCA-enrolled cases. Paired only with `peg1_controls`. |
| `peg_gwas` | In-house PEG GWAS (step 8). Paired only with `ppmi_controls`. |

---

## Step 1 — `01_exposure_candidates.R` → `exposure/`

### `exposure/{exposure}.tsv.gz` — candidate instruments
One row per (SNP, CpG) cis-mQTL pair with p < 5×10⁻⁸, after dropping SNP-in-probe pairs.

| Column | Description |
|---|---|
| `dataset` | Exposure name. |
| `cpg_id`, `cpg_chr`, `cpg_pos` | CpG and its position. |
| `snp_id`, `chr`, `pos` | SNP and its position. |
| `dist_bp` | SNP position − CpG position (pairs with \|dist\| ≤ 50 bp were removed as probe artefacts). |
| `ref`, `alt` | Alleles; `alt` is the effect allele. For meta exposures, both cohorts must have identical alleles. |
| `eaf` | ALT frequency (cohort `.afreq`; meta = sample-size-weighted mean). |
| `beta`, `se` | mQTL effect per ALT on the beta-value scale, and its SE. |
| `z` | Meta: meta z. Single cohort: MatrixEQTL t converted to z via its p-value (df = n − 2 − n_cov), so z is comparable across datasets. |
| `neglog10p` | −log10 p from z (computed without underflow). |
| `F_stat` | First-stage F for the single SNP = z². |
| `n` | Total samples in the exposure dataset. |
| `Q_pval`, `I2` | Cross-cohort heterogeneity of the mQTL effect (meta only; NA for single cohorts). |
| `flag_region` | See conventions. |

### `exposure/summary.tsv`
| Column | Description |
|---|---|
| `dataset`, `n` | Exposure and its sample size. |
| `pairs_p`, `cpgs_p` | Pairs / CpGs at p < 5×10⁻⁸ before filters. |
| `dropped_snp_in_probe` | Pairs with SNP ≤ 50 bp from the CpG. |
| `dropped_allele_mismatch` | SNPs whose REF/ALT differ between the meta's cohorts (multiallelic sites). |
| `dropped_no_allele_info` | Pairs whose SNP had no allele record. |
| `pairs_kept`, `cpgs_kept`, `snps_kept` | What remains. |
| `min_F` | Weakest F among kept pairs. |

---

## Step 2 — `02_harmonise_outcome.R` → `harmonised/`

### `harmonised/{exposure}__{outcome}.tsv.gz` — exposure SNPs matched to the GWAS
One row per exposure SNP that was found in the GWAS and kept.

| Column | Description |
|---|---|
| `snp_id`, `chr`, `pos`, `ref`, `alt` | SNP and exposure alleles. |
| `eaf_x` | Exposure ALT frequency. |
| `beta_y`, `se_y` | GWAS log-OR and SE, re-signed to the exposure ALT allele. |
| `eaf_y` | GWAS frequency of the exposure ALT allele. |
| `p_y` | GWAS p-value. |
| `harmonise` | How the GWAS alleles matched: `aligned` (effect allele = ALT), `swapped` (effect allele = REF; sign flipped), `strand_aligned` / `strand_swapped` (same after complementing). |
| `palindromic` | TRUE for A/T or C/G SNPs. Kept only if exposure MAF < 0.42, oriented by allele frequency. |

SNPs were dropped if palindromic and ambiguous, or if \|eaf_x − eaf_y\| > 0.2.

### `harmonised/summary.tsv`
| Column | Description |
|---|---|
| `exposure`, `outcome` | Combo. |
| `snps_exposure` | Unique exposure SNPs. |
| `at_gwas_position` | Exposure SNPs with any GWAS row at the same chr:pos. |
| `allele_matched` | Of those, alleles matched (directly or by strand). |
| `allele_mismatch` | At the position but different alleles. |
| `strand_flipped` | Matched only after complementing. |
| `palindromic_kept` | Palindromic SNPs kept (resolved by frequency). |
| `dropped_palindromic`, `dropped_eaf_discordant` | Removed by the two filters. |
| `harmonised` | Rows written. |
| `pct_missing_from_gwas` | 100 × (1 − harmonised / snps_exposure). |

---

## Step 3 — `03_clump.R` → `instruments/`

### `instruments/{exposure}__{outcome}.tsv.gz` — clumped instruments
One row per instrument (SNP) per CpG after greedy LD clumping **within each CpG** (r² < 0.001, 1000G EUR). Contains all step 1 columns (`dataset` … `flag_region`, see above) and all step 2 columns (`eaf_x` … `palindromic`), plus:

| Column | Description |
|---|---|
| `ld_ref` | Per CpG: `lead_in_ref` (the most significant SNP is in the LD reference and is the first instrument), `lead_not_in_ref` (top SNP missing from 1000G EUR; the best SNP that is present was used), `missing` (no candidate in the reference; top SNP kept unclumped). |
| `n_candidates` | Harmonised candidate SNPs for this CpG before clumping. |
| `n_candidates_in_ref` | Of those, present and polymorphic in 1000G EUR. |
| `n_instruments` | Instruments kept for this CpG. |

### `instruments/summary.tsv`
| Column | Description |
|---|---|
| `combo` | `{exposure}__{outcome}`. |
| `cpgs`, `instruments` | CpGs and total instrument rows. |
| `cpgs_1_iv`, `cpgs_2_iv`, `cpgs_3plus_iv` | CpGs by number of instruments. |
| `lead_in_ref`, `lead_not_in_ref`, `no_snp_in_ref` | CpGs by `ld_ref` category. |

`tmp_ld/` holds scratch files from clumping (SNP lists, PLINK logs). Not a result.

---

## Step 4 — `04_mr.R` → `mr/`

### `mr/{exposure}__{outcome}.tsv.gz` — MR estimate per CpG
| Column | Description |
|---|---|
| `exposure`, `outcome` | Combo. |
| `cpg_id`, `cpg_chr`, `cpg_pos` | CpG. |
| `method` | `wald` (1 instrument) or `ivw_mre` (≥ 2; IVW with multiplicative random effects). |
| `n_iv` | Instruments used. |
| `b`, `se`, `z`, `p` | MR estimate (log-OR per 1.0 beta-value), SE, z, two-sided p. Wald SE is first-order (se_y / \|bx\|). |
| `fdr`, `p_bonf` | BH FDR and Bonferroni p across CpGs in this analysis. |
| `or_per_0.1`, `or_per_0.1_lo`, `or_per_0.1_hi` | OR per 0.1 beta-value with 95% CI. |
| `se_2nd` | Second-order Wald SE (adds exposure uncertainty); Wald only. |
| `Q`, `Q_p` | Cochran's Q heterogeneity across instruments (≥ 2). |
| `b_wm`, `se_wm` | Weighted median (≥ 3 instruments). |
| `b_egger`, `se_egger`, `egger_int`, `egger_int_p` | MR-Egger slope, SE, intercept (pleiotropy test), intercept p (≥ 3). |
| `snps` | Instruments, `;`-separated. |
| `min_F`, `mean_F` | First-stage F across instruments. |
| `min_dist_bp` | Closest instrument–CpG distance. |
| `r2_x`, `r2_y` | Variance explained by the instruments in methylation (from z and n) and in PD liability. |
| `flag_region` | Region flag(s) of the instruments. |
| `ld_ref` | Clumping status (see step 3). |
| `lead_bx`, `lead_by`, `lead_p_y` | Exposure beta, GWAS beta and GWAS p of the strongest instrument. |
| `steiger_dir_ok` | TRUE if r2_x > r2_y (supports SNP → methylation → PD direction). |

### `mr/summary.tsv`
| Column | Description |
|---|---|
| `analysis` | Combo. |
| `cpgs`, `wald`, `ivw` | CpGs tested, by method. |
| `fdr05`, `bonf05` | CpGs at FDR < 0.05 / Bonferroni < 0.05. |
| `fdr05_flag_region` | FDR hits in a flagged region. |
| `fdr05_Q_het` | FDR hits with heterogeneity Q_p < 0.05. |
| `steiger_ok_pct` | % of CpGs with `steiger_dir_ok`. |
| `lambda_gc` | Genomic inflation of the per-CpG MR z². |

---

## Step 5 — `05_regional_mqtl.R` → `regional_mqtl/`

### `regional_mqtl/{cohort}.tsv.gz` — untruncated regional mQTL statistics
Every SNP within ±1 Mb of each target CpG (CpGs with FDR < 0.05 in any `*__nalls2019` MR), refitted with the PD-QTLs model. Cohorts: `peg1_controls`, `peg1_cases`, `ppmi_controls`, `ppmi_cases`.

| Column | Description |
|---|---|
| `cpg_id`, `snp_id` | Pair. |
| `dist_bp` | SNP − CpG position. |
| `beta`, `se`, `t`, `p` | Methylation beta-value per ALT allele, SE, t (df = n − 2 − n_cov), p. |
| `eaf` | ALT frequency in the analysed samples. |
| `cohort`, `n` | Cohort and sample size. |

### `regional_mqtl/validation.tsv` — refit vs PD-QTLs cis file on shared pairs
| Column | Description |
|---|---|
| `cohort`, `n` | Cohort. |
| `shared_pairs` | Pairs present in both. |
| `cor_beta`, `max_abs_diff_beta` | Agreement of beta. |
| `cor_t`, `max_rel_diff_t` | Agreement of t; the script stops if max_rel_diff_t ≥ 1e-3. |

---

## Step 6 — `06_coloc.R` → `coloc/`

### `coloc/{exposure}__{outcome}.tsv.gz` — coloc.abf per MR hit
One row per CpG with MR FDR < 0.05 in that combo. Exposure statistics: the regional mQTL (meta exposures re-meta-analysed per SNP over SNPs with identical alleles in both cohorts).

| Column | Description |
|---|---|
| `cpg_id` … `flag_region` | Copied from the MR file (`n_iv`, `b`, `se`, `p`, `fdr`, `or_per_0.1`). |
| `nsnps` | SNPs in the coloc region (regional mQTL ∩ harmonised GWAS). |
| `PP.H0` … `PP.H4` | Posteriors (p1 = p2 = 1e-4, p12 = 5e-6): H0 no signal; H1 mQTL only; H2 GWAS only; **H3 two distinct causal variants**; **H4 one shared causal variant**. |
| `PP.H4_p12_1e5` | PP.H4 with p12 = 1e-5 (sensitivity). |
| `top_snp_H4`, `top_snp_PP` | SNP with the highest posterior of being the shared variant, and that posterior. |
| `lead_mqtl_snp`, `lead_mqtl_z` | Strongest regional mQTL SNP and its z. |
| `lead_gwas_snp`, `lead_gwas_p` | Strongest GWAS SNP in the region and its p. |
| `gwas_p_at_lead_mqtl` | GWAS p at the lead mQTL SNP. |
| `coloc` | Call: `strong_H4` (H4 ≥ 0.8), `moderate_H4` (≥ 0.5), `distinct_H3` (H3 ≥ 0.5), `inconclusive`, `not_tested`. |
| `exposure`, `outcome` | Combo. |

### `coloc/summary.tsv`
| Column | Description |
|---|---|
| `exposure`, `outcome` | Combo. |
| `mr_hits`, `tested` | FDR hits and how many had coloc run. |
| `strong_H4`, `moderate_H4`, `distinct_H3`, `inconclusive` | Counts by call. |
| `median_nsnps` | Median region size. |

---

## Step 7 — `07_coloc_susie.R` → `coloc_susie/`

Files: `meta_controls__nalls2019` (primary), `meta_controls__nalls2019_noproxy` (no-proxy GWAS on its own MR hits), `meta_controls__nalls2019_noproxy_on_nalls2019` (no-proxy GWAS on the primary hits; sensitivity). mQTL side: PEG1 controls regional statistics with in-sample LD. GWAS side: 1000G EUR LD.

### `coloc_susie/{name}.tsv.gz` — one row per CpG
| Column | Description |
|---|---|
| `cpg_id`, `cpg_chr`, `cpg_pos`, `b`, `p`, `fdr`, `flag_region` | From the MR file that defined the hit set. |
| `nsnps` | SNPs used (regional ∩ GWAS ∩ 1000G EUR, polymorphic). |
| `cs_mqtl`, `cs_gwas` | Number of SuSiE 95% credible sets for each trait (NA if SuSiE failed). |
| `s_lambda_mqtl`, `s_lambda_gwas` | `estimate_s_rss` z-vs-LD consistency; near 0 is consistent, larger values indicate LD mismatch. |
| `best_H4`, `best_H3` | coloc.susie posteriors for the credible-set pair with the highest H4. |
| `hit_mqtl`, `hit_gwas` | Top SNPs of that pair's credible sets. |
| `abf_H3`, `abf_H4`, `abf_call` | coloc.abf result from step 6 for the same outcome (NA if the CpG was not a hit for that outcome). |
| `susie_call` | `strong_H4` / `moderate_H4` / `distinct_H3` / `inconclusive`; `no_credible_set` (a trait had no credible set, so coloc.susie can't test it); `not_run`. |

### `coloc_susie/{name}_pairs.tsv.gz` — every credible-set pair
| Column | Description |
|---|---|
| `nsnps` | SNPs. |
| `hit1`, `hit2` | Top SNP of the mQTL (1) and GWAS (2) credible set. |
| `PP.H0.abf` … `PP.H4.abf` | coloc posteriors for this pair. |
| `idx1`, `idx2` | Credible-set indices. |
| `cpg_id` | CpG. |

---

## Step 8 — `08_inhouse_gwas.R` → `inhouse_gwas/`

### `inhouse_gwas/{ppmi,ppmi_nogen,peg}.tsv.gz` — in-house PD GWAS (Nalls format)
PLINK2 logistic (Firth fallback), European ancestry, covariates sex, age, PC1–5; imputation R² ≥ 0.3, MAF ≥ 0.01.

| Column | Description |
|---|---|
| `chromosome`, `base_pair_location` | Variant position. |
| `effect_allele`, `other_allele` | Tested allele (PLINK A1) and the other. |
| `beta`, `standard_error` | log-OR per effect allele and SE. |
| `effect_allele_frequency` | A1 frequency in cases + controls. |
| `p_value` | p. |
| `n_cases`, `n_controls` | Sample counts (constant per file). |

**Individual-level, local only:** `*.pheno` and `*.covar` are the PLINK input files (per-participant phenotype and covariates). They must never be committed or shared. `*.log` / `*.plink.log` are PLINK logs.

---

## Step 9 — `09_results_table.R` → `tables/` (run after steps 10, 11 and 12)

### `tables/primary_hits.tsv` — one row per primary hit (meta_controls × nalls2019, FDR < 0.05)

**Primary MR and instrument**

| Column | Description |
|---|---|
| `cpg_id`, `cpg_chr`, `cpg_pos` | CpG. |
| `method`, `n_iv`, `b`, `se`, `p`, `fdr`, `or_per_0.1`, `or_per_0.1_lo`, `or_per_0.1_hi`, `steiger_dir_ok`, `flag_region` | From the primary MR file (step 4). |
| `iv_snp`, `iv_ref`, `iv_alt`, `iv_eaf` | Strongest instrument and its alleles / ALT frequency. |
| `iv_bx`, `iv_F`, `iv_dist_bp` | Its mQTL beta, F-statistic, and distance to the CpG. |

**Sensitivity and crossed MR** (`b_*` = MR estimate, `p_*` = p; NA if the CpG wasn't tested there)

| Suffix | Analysis |
|---|---|
| `_noproxy` | meta_controls × no-proxy GWAS. |
| `_casesmeta` | meta_cases × nalls2019. |
| `_peg1ctrl`, `_ppmictrl` | Single-cohort controls mQTL × nalls2019. |
| `_x_pegmqtl_ppmigwas`, `_x_pegmqtl_ppminogen` | Crossed: PEG1 controls mQTL × PPMI GWAS (with / without genetic-cohort cases). |
| `_x_ppmimqtl_peggwas` | Crossed: PPMI controls mQTL × PEG GWAS. |
| `crossed_tested` | Crossed analyses in which the CpG was tested (0–3). |
| `crossed_same_dir` | Of those, how many have the same sign as the primary `b`. |
| `crossed_nominal` | How many have the same sign and p < 0.05. |

**Colocalization**

| Column | Description |
|---|---|
| `abf_H3`, `abf_H4`, `abf_call` | coloc.abf (step 6). |
| `susie_cs_mqtl`, `susie_cs_gwas`, `susie_H4`, `susie_H3`, `susie_call` | coloc-SuSiE, primary GWAS (step 7). |
| `susie_noproxy_cs_gwas`, `susie_noproxy_H4`, `susie_noproxy_H3`, `susie_noproxy_call` | coloc-SuSiE with the no-proxy GWAS (sensitivity; does not affect the tier). |

**MAPT (step 10; NA/empty outside the MAPT region)**

| Column | Description |
|---|---|
| `mapt_class` | See step 10. |
| `mapt_r2_h2` | r² of the instrument with the H2 tag (1000G EUR). |
| `mapt_cond_mqtl_p` | p of the strongest mQTL left after conditioning on H2. |
| `mapt_cond_mr_p`, `mapt_cond_gwas_p` | H2-conditional MR p, and conditional GWAS p at that SNP. |
| `mapt_cond_abf_H4` | coloc.abf H4 on the H2-conditional statistics. |

**Shared instruments (step 11)**

| Column | Description |
|---|---|
| `iv_group`, `iv_group_size`, `iv_group_rank` | Instrument group, its size, and this CpG's rank in it (reporting order only). |
| `gene`, `gene_feature`, `island` | 450K manifest annotation (see step 11). |

**Probe-body SNPs (step 12; reported, does not change the tier)**

| Column | Description |
|---|---|
| `probe_type` | Infinium I or II. |
| `n_probe_snps`, `probe_snps` | SNPs inside the 50-mer probe body. |
| `probe_snp_max_r2_iv` | Highest EUR r² between a probe-body SNP and an instrument of the CpG. |
| `probe_snp_is_regional_lead` | TRUE if a probe-body SNP is the strongest regional mQTL SNP. |
| `probe_snp_ld` | TRUE if `probe_snp_max_r2_iv` ≥ 0.8: the mQTL may be partly a hybridisation artefact. |

**Tier and required metadata**

| Column | Description |
|---|---|
| `tier` | **A**: SuSiE H4 ≥ 0.8, or abf H4 ≥ 0.8 where SuSiE had no credible set. **B**: moderate (H4 0.5–0.8), same rule. **C**: everything else (distinct signals or no support). MAPT CpGs are C unless an H2-conditional coloc has H4 ≥ 0.8. |
| `exposure_sample`, `methylation_timepoint`, `mqtl_covariates`, `outcome_sample`, `estimand`, `instrument_rule`, `iv_assumptions` | Text, same in every row: the sample, timepoint, covariates, effect scale, instrument rule and the IV assumptions relied on (required by CLAUDE.md on every result table). |

### `tables/loci.tsv` — primary hits grouped into loci
A new locus starts when the next CpG (sorted by position) is on another chromosome or more than 1 Mb away.

| Column | Description |
|---|---|
| `locus` | Locus number. |
| `chr`, `start`, `end` | Span of the hit CpGs. |
| `n_cpg` | Hit CpGs in the locus. |
| `best_p`, `best_tier` | Smallest MR p and best tier. |
| `cpgs_tier_A`, `cpgs_tier_B` | CpGs by tier. |
| `flag_region` | Region flag. |
| `crossed_nominal_any` | Any CpG with `crossed_nominal` > 0. |

### `tables/instrument_groups.tsv`
| Column | Description |
|---|---|
| `iv_group` | Instrument group (step 11). |
| `n_cpg` | Primary hits in the group. |
| `best_tier` | Best tier among them. |
| `probe_snp_ld_cpgs` | CpGs in the group with `probe_snp_ld` (step 12). |

---

## Step 10 — `10_mapt.R` → `mapt/`

H2 tag: rs8070723, 17:44081064 A>G; **G (ALT) = H2**. Analyses use PEG1 + PPMI controls (meta) and the primary GWAS. Covers MAPT-region CpGs with FDR < 0.05 in any `*__nalls2019` MR.

### `mapt/cpgs.tsv`
| Column | Description |
|---|---|
| `cpg_id`, `cpg_pos` | CpG (chr 17). |
| `iv_snp`, `iv_r2_h2_eur` | Primary instrument and its r² with the H2 tag (NA if the CpG isn't a primary hit). |
| `marg_lead_snp`, `marg_lead_z`, `marg_lead_r2_h2` | Strongest marginal controls-meta mQTL SNP, its z, its r² with H2. |
| `h2_mqtl_beta`, `h2_mqtl_z` | mQTL effect of the H2 allele on this CpG. |
| `h2_wald_b`, `h2_wald_se` | Haplotype-level Wald ratio using the H2 tag as instrument. |
| `cond_lead_snp`, `cond_lead_z`, `cond_lead_p` | Strongest mQTL after adding H2 dosage as a covariate (exact, in-sample). |
| `cond_n_gw` | SNPs with conditional mQTL p < 5×10⁻⁸. |
| `cond_mr_snp`, `cond_mr_bx`, `cond_mr_px` | Conditional mQTL lead among SNPs available in the conditional GWAS, its beta and p. |
| `cond_mr_b`, `cond_mr_se`, `cond_mr_p` | H2-conditional Wald ratio at that SNP (NA if `cond_mr_px` ≥ 5×10⁻⁸). |
| `cond_gwas_p` | H2-conditional GWAS p at that SNP (approximate, 1000G EUR LD). |
| `cond_abf_H3`, `cond_abf_H4` | coloc.abf on conditional mQTL vs conditional GWAS (NA if no conditional mQTL at p < 5×10⁻⁸). |
| `r2_h2` | r² used for classification: `iv_r2_h2_eur`, or `marg_lead_r2_h2` when there's no primary instrument. |
| `mapt_class` | `H2_haplotype` (r² ≥ 0.5 and no conditional mQTL left), `H2_tagging_residual` (r² ≥ 0.5 but a conditional mQTL remains), `not_H2` (r² < 0.5), `lead_not_in_ref`. |
| `primary_b`, `primary_p`, `primary_fdr` | Primary MR result (NA if not a primary hit). |

### `mapt/gwas_h2.tsv` — H2 effect in each GWAS
| Column | Description |
|---|---|
| `outcome`, `n_cases`, `n_controls` | GWAS. |
| `h2_logOR`, `se`, `p` | log-OR per H2 (G) allele. |
| `h2_freq` | H2 frequency in the GWAS sample. |

### `mapt/breadth.tsv` — how many CpGs the H2 tag regulates
| Column | Description |
|---|---|
| `cohort` | mQTL cohort. |
| `cpgs_h2_p5e8`, `cpgs_h2_p1e5` | CpGs whose cis-mQTL includes the H2 tag at p < 5×10⁻⁸ / 1×10⁻⁵ (PD-QTLs cis files). |
| `n`, `h2_freq` | Sample size and H2 frequency (controls only). |

---

## Step 11 — `11_shared_instruments.R` → `shared_iv/`

Groups primary hits whose lead instruments are the same SNP or in LD r² ≥ 0.8 (1000G EUR, single linkage).

### `shared_iv/cpgs.tsv` — one row per primary hit
| Column | Description |
|---|---|
| `cpg_id`, `cpg_chr`, `cpg_pos` | CpG. |
| `group_id`, `group_size` | Instrument group (G01…, numbered in chromosome order (sorted as text, so 1, 10, 11, …), then position) and its size (1 = not shared). |
| `rank_in_group` | Reporting order: best PD coloc (`pd_H4`), then strongest mQTL. **Not** evidence of causality within the group. |
| `iv_snp`, `iv_z` | Lead instrument and its mQTL z. |
| `n_instruments_cpg` | Clumped instruments for this CpG. |
| `gene`, `gene_feature` | 450K manifest UCSC RefGene name(s) and feature(s) (TSS1500, TSS200, 5'UTR, 1stExon, Body, 3'UTR), de-duplicated; empty = intergenic. |
| `island` | Relation to CpG island: Island, N_/S_Shore, N_/S_Shelf, OpenSea. |
| `pd_H4` | SuSiE H4 if available, else coloc.abf H4. |

### `shared_iv/groups.tsv` — one row per multi-CpG group
| Column | Description |
|---|---|
| `group_id`, `chr`, `start`, `end` | Group and span. |
| `n_cpg` | CpGs. |
| `n_distinct_iv` | Distinct instrument SNPs across the group's CpGs. |
| `n_independent_iv` | Of those, mutually independent (r² < 0.01). |
| `mvmr_feasible` | TRUE if `n_independent_iv` ≥ `n_cpg` (needed for multivariable MR). |
| `genes` | Genes annotated to the group's CpGs. |
| `representative` | `rank_in_group` = 1 CpG (reporting only). |
| `n_pairs` | CpG pairs. |
| `min_mqtl_H4`, `median_mqtl_H4` | mQTL–mQTL coloc H4 across pairs (high = the CpGs share one mQTL causal variant). |
| `median_cor_peg1`, `median_cor_ppmi` | Median pairwise correlation of covariate-adjusted methylation in PEG1 / PPMI controls (co-methylation). |
| `flag_region` | Region flag. |

### `shared_iv/pairs.tsv` — one row per CpG pair within a group
| Column | Description |
|---|---|
| `group_id`, `cpg1`, `cpg2` | Pair. |
| `dist_bp` | Distance between the CpGs. |
| `same_iv` | TRUE if both have the identical lead instrument. |
| `cor_peg1`, `cor_ppmi` | Covariate-adjusted methylation correlation in each controls cohort. |
| `mqtl_H3`, `mqtl_H4` | mQTL–mQTL coloc.abf (controls meta): distinct vs shared mQTL causal variant. |

---

## Step 12 — `12_probe_snps.R` → `probe_snps/` (run before step 9)

Primary hits only. Each 450K/EPIC 50-mer (minfi manifest `ProbeSeqA`) is matched to the bisulfite-converted hg19 reference around its CpG (both strands, ≤ 3 mismatches, exactly one match required). Candidate SNPs are those in the step 5 regional files or in 1000G EUR. LD with the instruments uses 1000G EUR.

### `probe_snps/cpgs.tsv` — one row per primary hit
| Column | Description |
|---|---|
| `cpg_id` | CpG. |
| `array` | Manifest used (450K, or EPIC for EPIC-only probes). |
| `probe_type` | Infinium I or II. |
| `probe_start`, `probe_end` | Probe body on hg19. |
| `sbe_pos` | Type I single-base-extension position (NA for Type II, whose extension base is the CpG itself). |
| `n_probe_snps`, `probe_snps` | SNPs inside the probe body. |
| `probe_snp_positions` | Where they sit: `body`, `CpG`, or `SBE` (the last two distort the beta-value most). |
| `probe_snp_max_r2_iv` | Highest EUR r² between a probe-body SNP and an instrument. |
| `probe_snp_min_mqtl_p` | Strongest regional mQTL p of a probe-body SNP (primary-exposure cohorts). |
| `probe_snp_is_regional_lead` | A probe-body SNP is the strongest regional mQTL SNP. |
| `probe_snp_ld` | `probe_snp_max_r2_iv` ≥ 0.8. |

### `probe_snps/snps.tsv` — one row per CpG × probe-body SNP
| Column | Description |
|---|---|
| `cpg_id`, `snp_id` | Pair. |
| `dist_bp`, `position`, `probe_type` | Distance to the CpG, position class, probe type. |
| `max_r2_iv`, `iv_best` | Highest r² with an instrument, and which instrument. |
| `in_1000g` | SNP is in the 1000G EUR reference (r² is NA otherwise). |
| `mqtl_p_min` | Strongest regional mQTL p (primary-exposure cohorts). |
| `is_regional_lead` | The SNP is the regional mQTL lead. |
