#!/usr/bin/env Rscript
# Phase 1, step 1: candidate cis-mQTL instruments per exposure dataset.
#
# For every (SNP, CpG) pair at p < instruments$p_threshold:
#   - attach CpG position and drop SNP-in-probe artefacts (|SNP - CpG| <= snp_in_probe_bp)
#   - attach REF/ALT from each pipeline's positions/export_allele.txt. mQTL beta is per
#     ALT (the counted allele), so effect_allele = ALT, other_allele = REF.
#   - meta datasets: require identical REF/ALT in all contributing pipelines (each
#     pipeline picks its multiallelic record by its own frequencies; mismatches are
#     a handful of multiallelic sites and are dropped)
#   - attach ALT frequency (cohort-wide .afreq; meta = n-weighted mean)
# Output: results/phase1/exposure/{dataset}.tsv.gz  (all candidates; clumping is step 3)
#
# Usage: Rscript scripts/R/01_exposure_candidates.R [dataset ...]   (default: all in config)

source("scripts/R/phase1_lib.R")
cfg  <- load_config()
args <- commandArgs(trailingOnly = TRUE)
dsets <- if (length(args)) args else names(cfg$exposures)
pthr  <- cfg$instruments$p_threshold
zthr  <- z_for_p(pthr)

# ── Read candidate pairs ─────────────────────────────────────────────────────
read_pairs <- function(e) {
  f <- pdqtl_path(cfg, e$file)
  if (e$type == "meta") {
    # pair_id,snp_id,cpg_id,meta_beta,meta_se,meta_pval,meta_z,...
    cmd <- sprintf("zcat %s | awk -F, 'NR==1 || $7 > %.6f || $7 < -%.6f' | cut -d, -f2-5,7-10",
                   shQuote(f), zthr, zthr)
    d <- fread(cmd = cmd, colClasses = list(character = c("snp_id", "cpg_id")))
    setnames(d, c("meta_beta", "meta_se", "meta_z"), c("beta", "se", "z"))
    d[, Q := NULL]
  } else {
    # MatrixEQTL: SNP gene beta t-stat p-value FDR ; p-sorted
    cmd <- sprintf("zcat %s | awk -F'\\t' 'NR==1 || $5 < %g' | cut -f1-4", shQuote(f), pthr)
    d <- fread(cmd = cmd, colClasses = list(character = c("SNP", "gene")))
    setnames(d, c("SNP", "gene", "t-stat"), c("snp_id", "cpg_id", "tstat"))
    d[, se := abs(beta / tstat)]
    d[, `:=`(Q_pval = NA_real_, I2 = NA_real_)]
  }
  d[, snp_id := sub(":[ACGT]+:[ACGT]+$", "", snp_id)]
  d[, c("chr", "pos") := tstrsplit(snp_id, ":", fixed = TRUE)]
  d[, pos := as.integer(pos)]
  d
}

# ── Allele + frequency tables (one per pipeline, cached in memory) ───────────
allele_cache <- list()
pipeline_alleles <- function(pipe, snps) {
  if (is.null(allele_cache[[pipe]])) {
    logmsg("  reading %s export_allele + afreq", pipe)
    ea <- fread(pdqtl_path(cfg, pipe, "positions/export_allele.txt"), header = FALSE,
                col.names = c("id", "counted"))
    af <- fread(pdqtl_path(cfg, pipe, "positions/cohort_allele_freq.afreq"),
                select = c("ID", "REF", "ALT", "ALT_FREQS"))
    x <- merge(ea, af, by.x = "id", by.y = "ID")
    stopifnot(nrow(x) == nrow(ea), all(x$counted == x$ALT))
    x[, snp_id := sub(":[ACGT]+:[ACGT]+$", "", id)]
    stopifnot(!anyDuplicated(x$snp_id))
    allele_cache[[pipe]] <<- x[, .(snp_id, ref = REF, alt = ALT, eaf = ALT_FREQS)]
  }
  allele_cache[[pipe]][snp_id %in% snps]
}

probe_pos <- function(e) {
  ps <- unique(vapply(e$cohorts, function(c) cfg$cohorts[[c]]$probes, ""))
  p <- rbindlist(lapply(ps, function(f) fread(pdqtl_path(cfg, f), colClasses = "character")))
  p <- unique(p)
  stopifnot(!anyDuplicated(p$cpg_id))
  p[, .(cpg_id, cpg_chr = chr, cpg_pos = as.integer(pos))]
}

# ── Main ─────────────────────────────────────────────────────────────────────
summary_rows <- list()
for (ds in dsets) {
  e <- cfg$exposures[[ds]]
  logmsg("=== %s (%s) ===", ds, e$type)
  d <- read_pairs(e)
  n_pairs0 <- nrow(d); n_cpg0 <- uniqueN(d$cpg_id)
  logmsg("  pairs at p < %g: %s (%s CpGs)", pthr, format(n_pairs0, big.mark = ","), format(n_cpg0, big.mark = ","))

  # For single cohorts MatrixEQTL gives t, not z. With n as small as 70 and
  # 11-13 covariates, t and z differ; convert t -> p (t on n - 2 - n_cov df)
  # -> z so that neglog10p and F are comparable across datasets.
  ns <- vapply(e$cohorts, function(c) cfg$cohorts[[c]]$n_samples, 0)
  n_total <- sum(ns)
  if (e$type == "cohort") {
    df <- n_total - 2 - cfg$cohorts[[e$cohorts[[1]]]]$n_cov
    log_p1 <- pt(-abs(d$tstat), df = df, log.p = TRUE)          # one-sided
    d[, z := sign(tstat) * qnorm(log_p1, lower.tail = FALSE, log.p = TRUE)]
    d[, tstat := NULL]
  }

  # CpG position + SNP-in-probe filter
  d <- merge(d, probe_pos(e), by = "cpg_id", all.x = TRUE)
  stopifnot(!anyNA(d$cpg_pos))
  d[, dist_bp := pos - cpg_pos]
  n_inprobe <- d[abs(dist_bp) <= cfg$instruments$snp_in_probe_bp, .N]
  d <- d[abs(dist_bp) > cfg$instruments$snp_in_probe_bp]

  # Alleles: all contributing pipelines must agree
  pipes <- unique(vapply(e$cohorts, function(c) cfg$cohorts[[c]]$pipeline, ""))
  snps <- unique(d$snp_id)
  al <- lapply(pipes, pipeline_alleles, snps = snps)
  names(al) <- pipes
  a <- al[[1]]
  n_allele_mismatch <- 0L
  if (length(pipes) > 1) {
    for (p in pipes[-1]) {
      b <- al[[p]]
      a <- merge(a, b[, .(snp_id, ref_b = ref, alt_b = alt, eaf_b = eaf)], by = "snp_id")
      bad <- a$ref != a$ref_b | a$alt != a$alt_b
      n_allele_mismatch <- n_allele_mismatch + sum(bad)
      a <- a[!bad]
      # n-weighted ALT frequency across pipelines
      n_a <- sum(ns[vapply(e$cohorts, function(c) cfg$cohorts[[c]]$pipeline, "") != p])
      n_b <- sum(ns[vapply(e$cohorts, function(c) cfg$cohorts[[c]]$pipeline, "") == p])
      a[, eaf := (eaf * n_a + eaf_b * n_b) / (n_a + n_b)]
      a[, c("ref_b", "alt_b", "eaf_b") := NULL]
    }
  }
  n_no_allele <- sum(!d$snp_id %in% a$snp_id)
  d <- merge(d, a, by = "snp_id")

  d[, `:=`(neglog10p = neglog10p_z(z), F_stat = z^2, n = n_total,
           flag_region = flag_region(chr, pos), dataset = ds)]
  setcolorder(d, c("dataset", "cpg_id", "cpg_chr", "cpg_pos", "snp_id", "chr", "pos", "dist_bp",
                   "ref", "alt", "eaf", "beta", "se", "z", "neglog10p", "F_stat", "n",
                   "Q_pval", "I2", "flag_region"))
  setorder(d, cpg_id, -neglog10p)
  fwrite(d, out_path(cfg, "exposure", paste0(ds, ".tsv.gz")), sep = "\t")

  summary_rows[[ds]] <- data.table(
    dataset = ds, n = n_total, pairs_p = n_pairs0, cpgs_p = n_cpg0,
    dropped_snp_in_probe = n_inprobe, dropped_allele_mismatch = n_allele_mismatch,
    dropped_no_allele_info = n_no_allele, pairs_kept = nrow(d), cpgs_kept = uniqueN(d$cpg_id),
    snps_kept = uniqueN(d$snp_id), min_F = round(min(d$F_stat), 1))
  logmsg("  kept %s pairs, %s CpGs, %s unique SNPs (in-probe dropped %s, allele mismatch %s, no allele info %s)",
         format(nrow(d), big.mark = ","), format(uniqueN(d$cpg_id), big.mark = ","),
         format(uniqueN(d$snp_id), big.mark = ","), n_inprobe, n_allele_mismatch, n_no_allele)
}
s <- rbindlist(summary_rows)
sf <- out_path(cfg, "exposure", "summary.tsv")
if (file.exists(sf) && length(args)) s <- rbind(fread(sf)[!dataset %in% s$dataset], s)
fwrite(s, sf, sep = "\t")
print(s)
