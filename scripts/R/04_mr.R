#!/usr/bin/env Rscript
# Phase 1, step 4: two-sample MR per CpG (exposure = methylation beta-value, outcome = PD log-OR).
#
# Estimand: log-OR of PD per 1-unit increase in methylation beta-value (0 -> 1). Because
# a full 0->1 change is not biologically meaningful, OR per 0.1 (10 percentage points)
# is also reported. Binary outcome on the log-OR scale (ratio methods); ORs are
# non-collapsible and this is a population-averaged two-sample estimate.
#
# Per CpG with k clumped instruments:
#   k = 1  Wald ratio  b = by/bx ; se = se_y/|bx| (first order; second-order se also given)
#   k >= 2 IVW, multiplicative random effects (se inflated by sqrt(Q/(k-1)) when > 1),
#          Cochran's Q for heterogeneity
#   k >= 3 + weighted median and MR-Egger (MendelianRandomization)
# Steiger directionality: sum r^2 of instruments with methylation (from z and n) vs with PD
# (liability approximation r^2 = vg/(vg + pi^2/3), vg = by^2 * 2p(1-p)). With strong cis
# instruments this essentially always favours SNP -> methylation -> PD; reported anyway.
# FDR: Benjamini-Hochberg across CpGs within each exposure x outcome analysis.
#
# Output: results/phase1/mr/{exposure}__{outcome}.tsv.gz, results/phase1/mr/summary.tsv
# Usage:  Rscript scripts/R/04_mr.R

source("scripts/R/phase1_lib.R")
suppressPackageStartupMessages(library(MendelianRandomization))
cfg <- load_config()

ifiles <- list.files(out_path(cfg, "instruments"), pattern = "__.*\\.tsv\\.gz$", full.names = TRUE)

mr_cpg <- function(bx, sx, by, sy) {
  k <- length(bx)
  if (k == 1) {
    b <- by / bx
    se <- sy / abs(bx)
    se2 <- sqrt(sy^2 / bx^2 + by^2 * sx^2 / bx^4)
    return(list(method = "wald", b = b, se = se, se_2nd = se2, Q = NA_real_, Q_p = NA_real_,
                b_wm = NA_real_, se_wm = NA_real_, b_egger = NA_real_, se_egger = NA_real_,
                egger_int = NA_real_, egger_int_p = NA_real_))
  }
  w  <- bx^2 / sy^2
  b  <- sum(bx * by / sy^2) / sum(w)
  se_fe <- 1 / sqrt(sum(w))
  Q  <- sum((by / bx - b)^2 * w)
  phi <- max(1, Q / (k - 1))
  out <- list(method = "ivw_mre", b = b, se = se_fe * sqrt(phi), se_2nd = NA_real_, Q = Q,
              Q_p = pchisq(Q, k - 1, lower.tail = FALSE),
              b_wm = NA_real_, se_wm = NA_real_, b_egger = NA_real_, se_egger = NA_real_,
              egger_int = NA_real_, egger_int_p = NA_real_)
  if (k >= 3) {
    inp <- mr_input(bx = bx, bxse = sx, by = by, byse = sy)
    wm <- tryCatch(mr_median(inp, weighting = "weighted", iterations = 1000), error = function(e) NULL)
    eg <- tryCatch(mr_egger(inp), error = function(e) NULL)
    if (!is.null(wm)) { out$b_wm <- wm@Estimate; out$se_wm <- wm@StdError }
    if (!is.null(eg)) { out$b_egger <- eg@Estimate; out$se_egger <- eg@StdError.Est
                        out$egger_int <- eg@Intercept; out$egger_int_p <- eg@Pvalue.Int }
  }
  out
}

summ <- list()
set.seed(20261115)
for (f in ifiles) {
  cb <- sub("\\.tsv\\.gz$", "", basename(f))
  oc <- sub(".*__", "", cb)
  d  <- fread(f, colClasses = list(character = c("snp_id", "chr", "cpg_chr")))
  logmsg("%s: %s CpGs, %s instruments", cb, format(uniqueN(d$cpg_id), big.mark = ","),
         format(nrow(d), big.mark = ","))

  d[, r2_x := z^2 / (z^2 + n - 2)]
  d[, r2_y := { vg <- beta_y^2 * 2 * eaf_y * (1 - eaf_y); vg / (vg + pi^2 / 3) }]

  r <- d[, {
    m <- mr_cpg(beta, se, beta_y, se_y)
    c(m, list(n_iv = .N, snps = paste(snp_id, collapse = ";"),
              min_F = min(F_stat), mean_F = mean(F_stat),
              min_dist_bp = min(abs(dist_bp)),
              r2_x = sum(r2_x), r2_y = sum(r2_y),
              flag_region = paste(unique(na.omit(flag_region)), collapse = ";"),
              ld_ref = ld_ref[1], cpg_chr = cpg_chr[1], cpg_pos = cpg_pos[1],
              lead_bx = beta[1], lead_by = beta_y[1], lead_p_y = p_y[1]))
  }, by = cpg_id]

  r[, z := b / se]
  r[, `:=`(p = 2 * pnorm(-abs(z)), or_per_0.1 = exp(0.1 * b),
           or_per_0.1_lo = exp(0.1 * (b - 1.96 * se)), or_per_0.1_hi = exp(0.1 * (b + 1.96 * se)),
           steiger_dir_ok = r2_x > r2_y)]
  r[, `:=`(fdr = p.adjust(p, "BH"), p_bonf = pmin(1, p * .N))]
  r[flag_region == "", flag_region := NA_character_]
  setorder(r, p)
  r[, `:=`(exposure = sub("__.*", "", cb), outcome = oc)]
  setcolorder(r, c("exposure", "outcome", "cpg_id", "cpg_chr", "cpg_pos", "method", "n_iv",
                   "b", "se", "z", "p", "fdr", "p_bonf", "or_per_0.1", "or_per_0.1_lo", "or_per_0.1_hi"))
  fwrite(r, out_path(cfg, "mr", paste0(cb, ".tsv.gz")), sep = "\t")

  summ[[cb]] <- data.table(
    analysis = cb, cpgs = nrow(r), wald = r[method == "wald", .N], ivw = r[method != "wald", .N],
    fdr05 = r[fdr < cfg$fdr, .N], bonf05 = r[p_bonf < 0.05, .N],
    fdr05_flag_region = r[fdr < cfg$fdr & !is.na(flag_region), .N],
    fdr05_Q_het = r[fdr < cfg$fdr & Q_p < 0.05, .N],
    steiger_ok_pct = round(100 * mean(r$steiger_dir_ok), 1),
    lambda_gc = round(median(r$z^2) / qchisq(0.5, 1), 3))
  logmsg("  FDR<%.2f: %d CpGs (Bonferroni: %d); lambda %.3f", cfg$fdr, summ[[cb]]$fdr05,
         summ[[cb]]$bonf05, summ[[cb]]$lambda_gc)
}
s <- rbindlist(summ)
fwrite(s, out_path(cfg, "mr", "summary.tsv"), sep = "\t")
print(s)
