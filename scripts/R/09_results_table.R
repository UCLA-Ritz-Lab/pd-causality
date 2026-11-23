#!/usr/bin/env Rscript
# Phase 1, step 9: consolidated results table for the primary MR hits.
#
# One row per primary-hit CpG (meta_controls x nalls2019, FDR < cfg$fdr) with:
#   primary MR estimate; sensitivity MR (no-proxy GWAS, cases-meta, single-cohort mQTLs);
#   crossed replication with in-house GWAS (PEG1-controls mQTL x PPMI GWAS, with and
#   without genetic-cohort cases; PPMI-controls mQTL x PEG GWAS); coloc.abf and
#   coloc-SuSiE; an evidence tier; and the metadata CLAUDE.md requires on every result
#   table (sample, cohort, methylation timepoint, covariates, instrument strength, IV
#   assumptions relied on).
# Summary-level only. Output: results/phase1/tables/primary_hits.tsv, loci.tsv
# Usage:  Rscript scripts/R/09_results_table.R

source("scripts/R/phase1_lib.R")
cfg <- load_config()
rd <- function(...) { f <- out_path(cfg, ...); if (file.exists(f)) suppressWarnings(fread(f, colClasses = list(character = "cpg_chr"))) else NULL }

p <- rd("mr", "meta_controls__nalls2019.tsv.gz")[fdr < cfg$fdr]
inst <- fread(out_path(cfg, "instruments", "meta_controls__nalls2019.tsv.gz"),
              colClasses = list(character = c("snp_id", "chr", "cpg_chr")))[cpg_id %in% p$cpg_id]
lead_iv <- inst[order(-neglog10p)][!duplicated(cpg_id), .(cpg_id, iv_snp = snp_id, iv_ref = ref, iv_alt = alt, iv_eaf = round(eaf, 3),
                                                         iv_bx = beta, iv_F = round(F_stat, 1), iv_dist_bp = dist_bp)]
out <- merge(p[, .(cpg_id, cpg_chr, cpg_pos, method, n_iv, b, se, p, fdr, or_per_0.1, or_per_0.1_lo, or_per_0.1_hi,
                   steiger_dir_ok, flag_region)], lead_iv, by = "cpg_id")

add_mr <- function(out, cb, tag) {
  x <- rd("mr", paste0(cb, ".tsv.gz"))
  if (is.null(x)) return(out)
  x <- x[, .(cpg_id, b, p)]
  setnames(x, c("b", "p"), paste0(c("b_", "p_"), tag))
  merge(out, x, by = "cpg_id", all.x = TRUE)
}
out <- add_mr(out, "meta_controls__nalls2019_noproxy", "noproxy")
out <- add_mr(out, "meta_cases__nalls2019", "casesmeta")
out <- add_mr(out, "peg1_controls__nalls2019", "peg1ctrl")
out <- add_mr(out, "ppmi_controls__nalls2019", "ppmictrl")
out <- add_mr(out, "peg1_controls__ppmi_gwas", "x_pegmqtl_ppmigwas")
out <- add_mr(out, "peg1_controls__ppmi_nogen_gwas", "x_pegmqtl_ppminogen")
out <- add_mr(out, "ppmi_controls__peg_gwas", "x_ppmimqtl_peggwas")

ab <- rd("coloc", "meta_controls__nalls2019.tsv.gz")[, .(cpg_id, abf_H3 = round(PP.H3, 3), abf_H4 = round(PP.H4, 3), abf_call = coloc)]
su <- rd("coloc_susie", "meta_controls__nalls2019.tsv.gz")[, .(cpg_id, susie_cs_mqtl = cs_mqtl, susie_cs_gwas = cs_gwas,
                                                              susie_H4 = round(best_H4, 3), susie_H3 = round(best_H3, 3), susie_call)]
out <- merge(merge(out, ab, by = "cpg_id", all.x = TRUE), su, by = "cpg_id", all.x = TRUE)
# sensitivity (does not change the tier): coloc-SuSiE with the no-proxy GWAS on the same hits
sn <- rd("coloc_susie", "meta_controls__nalls2019_noproxy_on_nalls2019.tsv.gz")
if (!is.null(sn)) out <- merge(out, sn[, .(cpg_id, susie_noproxy_cs_gwas = cs_gwas, susie_noproxy_H4 = round(best_H4, 3),
                                           susie_noproxy_H3 = round(best_H3, 3), susie_noproxy_call = susie_call)],
                               by = "cpg_id", all.x = TRUE)

# MAPT H1/H2 (step 10): H2-tag LD of the instrument, and H2-conditional MR / coloc
mp <- rd("mapt", "cpgs.tsv")
if (!is.null(mp)) out <- merge(out, mp[, .(cpg_id, mapt_class, mapt_r2_h2 = round(r2_h2, 3), mapt_cond_mqtl_p = cond_lead_p,
                                           mapt_cond_mr_p = cond_mr_p, mapt_cond_gwas_p = cond_gwas_p,
                                           mapt_cond_abf_H4 = round(cond_abf_H4, 3))], by = "cpg_id", all.x = TRUE)

# Shared instruments (step 11): CpGs in one group are one nomination, not separate ones
sg <- rd("shared_iv", "cpgs.tsv")
if (!is.null(sg)) out <- merge(out, sg[, .(cpg_id, iv_group = group_id, iv_group_size = group_size, iv_group_rank = rank_in_group,
                                           gene, gene_feature, island)], by = "cpg_id", all.x = TRUE)

# Probe-body SNPs (step 12): reported next to the tier, does not change it
ps <- rd("probe_snps", "cpgs.tsv")
if (!is.null(ps)) out <- merge(out, ps[, .(cpg_id, probe_type, n_probe_snps, probe_snps,
                                           probe_snp_max_r2_iv = round(probe_snp_max_r2_iv, 3),
                                           probe_snp_is_regional_lead, probe_snp_ld)], by = "cpg_id", all.x = TRUE)

# crossed replication: same direction as primary and p < 0.05 in any crossed analysis
xcols <- grep("^b_x_", names(out), value = TRUE)
out[, crossed_same_dir := rowSums(sapply(xcols, function(cn) sign(get(cn)) == sign(b)), na.rm = TRUE)]
out[, crossed_tested := rowSums(!is.na(.SD)), .SDcols = xcols]
out[, crossed_nominal := rowSums(sapply(xcols, function(cn) sign(get(cn)) == sign(b) & get(sub("^b_", "p_", cn)) < 0.05), na.rm = TRUE)]

# Evidence tier (coloc is the exclusion-restriction / LD-confounding check):
#   A: shared variant supported by coloc-SuSiE (H4 >= 0.8), or by coloc.abf (>= 0.8) where
#      SuSiE could not run, and no SuSiE H3 >= 0.5
#   B: moderate shared-variant support (H4 0.5-0.8 by either method, no SuSiE H3 >= 0.5)
#   C: distinct signals (H3 >= 0.5 by SuSiE, or by abf where SuSiE had no credible set)
out[, tier := fcase(
  susie_call == "strong_H4", "A",
  susie_call %in% c("no_credible_set", "not_run") & abf_call == "strong_H4", "A",
  susie_call == "moderate_H4", "B",
  susie_call %in% c("no_credible_set", "not_run") & abf_call == "moderate_H4", "B",
  default = "C")]
# MAPT: no CpG-attributable estimate unless an H2-conditional signal colocalizes
if ("mapt_class" %in% names(out))
  out[!is.na(mapt_class) & (is.na(mapt_cond_abf_H4) | mapt_cond_abf_H4 < 0.8), tier := "C"]

# metadata required on every result table
out[, `:=`(
  exposure_sample = "PEG1 controls (n=231, 450K) + PPMI controls (n=69, EPIC v1); common-effect IVW meta of per-cohort cis-mQTLs",
  methylation_timepoint = "baseline blood (controls: no PD diagnosis, no levodopa)",
  mqtl_covariates = "PEG: sex, age, 6 cell types, 5 ancestry PCs; PPMI: age, 5 cell types, 5 ancestry PCs (no sex)",
  outcome_sample = "Nalls 2019 GCST009325: 15,056 PD + 18,618 UKB proxy cases vs 449,056 controls, European, excl. 23andMe",
  estimand = "log-OR of PD per 1.0 beta-value (OR per 0.1 given); two-sample, ratio method",
  instrument_rule = "cis (<=1 Mb) p<5e-8, SNP >50 bp from CpG, clumped per CpG r2<0.001 (1000G EUR)",
  iv_assumptions = paste("relevance: F reported;",
                         "independence: two-sample, ancestry-matched (EUR GWAS; mQTL adjusted for 5 PCs); small PPMI-in-Nalls overlap;",
                         "exclusion: single instrument for most CpGs, pleiotropy-robust MR not possible -> relies on coloc;",
                         "neighbouring CpGs sharing one instrument are not distinguishable;",
                         "measurement: probe_snp_ld = instrument in LD (r2>=0.8) with a SNP in the probe body, so the mQTL may be partly a hybridisation artefact"))]
setorder(out, tier, p)
fwrite(out, out_path(cfg, "tables", "primary_hits.tsv"), sep = "\t")

setorder(out, cpg_chr, cpg_pos)
out[, locus := cumsum(c(TRUE, diff(cpg_pos) > 1e6 | cpg_chr[-1] != cpg_chr[-.N]))]
loci <- out[, .(chr = cpg_chr[1], start = min(cpg_pos), end = max(cpg_pos), n_cpg = .N,
                best_p = min(p), best_tier = min(tier), cpgs_tier_A = sum(tier == "A"), cpgs_tier_B = sum(tier == "B"),
                flag_region = flag_region[1], crossed_nominal_any = any(crossed_nominal > 0)), by = locus]
setorder(loci, best_tier, best_p)
fwrite(loci, out_path(cfg, "tables", "loci.tsv"), sep = "\t")
print(out[, .N, by = tier][order(tier)])
if ("iv_group" %in% names(out)) {
  ig <- out[, .(n_cpg = .N, best_tier = min(tier),
                probe_snp_ld_cpgs = if ("probe_snp_ld" %in% names(out)) sum(probe_snp_ld, na.rm = TRUE) else NA_integer_), by = iv_group]
  fwrite(ig, out_path(cfg, "tables", "instrument_groups.tsv"), sep = "\t")
  print(ig[, .(groups = .N, cpgs = sum(n_cpg)), by = best_tier][order(best_tier)])
}
print(loci)
