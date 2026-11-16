#!/usr/bin/env Rscript
# Phase 1, step 6: colocalization of the methylation (mQTL) and PD (GWAS) signals at each
# MR-significant CpG.
#
# A single-SNP Wald ratio cannot tell a shared causal variant (H4) from two distinct
# causal variants in LD (H3), which is the main LD-confounding threat to these MR
# estimates. coloc.abf (single causal variant per trait) is run on the untruncated
# regional statistics from step 5 (all SNPs within +-1 Mb of the CpG):
#   exposure: the same mQTL dataset as the MR analysis. Meta datasets are re-meta-analysed
#             per SNP (common-effect IVW) over SNPs present, with identical alleles, in
#             every contributing cohort; ALT frequency is the n-weighted mean.
#   outcome:  Nalls GWAS harmonised to the ALT allele (same rules as step 2).
# Priors: p1 = p2 = 1e-4; p12 = 5e-6 (primary, Wallace 2020) and 1e-5 (coloc default,
# reported as sensitivity). Loci with substantial H3, or several independent signals,
# are candidates for coloc-SuSiE (not run here).
#
# Output: results/phase1/coloc/{exposure}__{outcome}.tsv.gz (one row per CpG), summary.tsv
# Usage:  Rscript scripts/R/06_coloc.R

source("scripts/R/phase1_lib.R")
suppressPackageStartupMessages(library(coloc))
cfg <- load_config()
ic  <- cfg$instruments

# ── Regional mQTL stats per cohort + alleles ─────────────────────────────────
reg <- list()
for (co in names(cfg$cohorts)) {
  f <- out_path(cfg, "regional_mqtl", paste0(co, ".tsv.gz"))
  if (file.exists(f)) reg[[co]] <- fread(f, colClasses = list(character = c("cpg_id", "snp_id")))
}
logmsg("regional mQTL cohorts available: %s", paste(names(reg), collapse = ", "))
reg_snps <- unique(unlist(lapply(reg, `[[`, "snp_id")))

alleles <- list()
for (pipe in unique(vapply(cfg$cohorts, `[[`, "", "pipeline"))) {
  ea <- fread(pdqtl_path(cfg, pipe, "positions/export_allele.txt"), header = FALSE, col.names = c("id", "counted"))
  ea[, snp_id := sub(":[ACGT]+:[ACGT]+$", "", id)]
  ea <- ea[snp_id %in% reg_snps]
  ea[, c("chr", "pos", "ref", "alt") := tstrsplit(id, ":", fixed = TRUE)]
  stopifnot(all(ea$alt == ea$counted))
  alleles[[pipe]] <- ea[, .(snp_id, chr, pos = as.integer(pos), ref, alt)]
}

exposure_regional <- function(ds) {
  e <- cfg$exposures[[ds]]
  cos <- unlist(e$cohorts)
  if (!all(cos %in% names(reg))) return(NULL)
  parts <- lapply(cos, function(co) {
    a <- alleles[[cfg$cohorts[[co]]$pipeline]]
    merge(reg[[co]][, .(cpg_id, snp_id, beta, se, eaf, n)], a, by = "snp_id")[, cohort := co]
  })
  if (length(parts) == 1) return(parts[[1]][, .(cpg_id, snp_id, chr, pos, ref, alt, beta, se, eaf, n)])
  d <- rbindlist(parts)
  # keep SNPs present in every cohort with identical alleles
  d[, ok := .N == length(cos) & uniqueN(paste(ref, alt)) == 1, by = .(cpg_id, snp_id)]
  d <- d[ok == TRUE]
  d[, w := 1 / se^2]
  d[, .(chr = chr[1], pos = pos[1], ref = ref[1], alt = alt[1],
        beta = sum(w * beta) / sum(w), se = sqrt(1 / sum(w)),
        eaf = sum(n * eaf) / sum(n), n = sum(n)), by = .(cpg_id, snp_id)]
}

run_coloc <- function(d, o, p12) {
  s <- o$n_cases / (o$n_cases + o$n_controls)
  maf <- pmin(d$eaf, 1 - d$eaf)
  capture.output(r <- suppressWarnings(suppressMessages(coloc.abf(
    dataset1 = list(beta = d$beta, varbeta = d$se^2, snp = d$snp_id, type = "quant",
                    N = d$n[1], MAF = maf),
    dataset2 = list(beta = d$beta_y, varbeta = d$se_y^2, snp = d$snp_id, type = "cc",
                    s = s, N = o$n_cases + o$n_controls, MAF = maf),
    p1 = 1e-4, p2 = 1e-4, p12 = p12))))
  list(pp = r$summary, top = r$results$snp[which.max(r$results$SNP.PP.H4)],
       top_pp = max(r$results$SNP.PP.H4))
}

summ <- list()
for (oc in names(cfg$outcomes)) {
  o <- cfg$outcomes[[oc]]
  logmsg("=== outcome %s", oc)
  g <- read_gwas(o$file)
  for (ds in names(cfg$exposures)) {
    mrf <- out_path(cfg, "mr", sprintf("%s__%s.tsv.gz", ds, oc))
    if (!file.exists(mrf)) next
    hits <- fread(mrf, colClasses = list(character = "cpg_chr"))[fdr < cfg$fdr]
    x <- exposure_regional(ds)
    if (is.null(x) || !nrow(hits)) next
    x <- x[cpg_id %in% hits$cpg_id]
    xs <- unique(x[, .(snp_id, chr, pos, ref, alt, eaf)])
    h <- harmonise_gwas(xs, g[unique(xs[, .(chr, pos)]), on = c("chr", "pos"), nomatch = NULL], ic)[is.na(drop)]
    x <- merge(x, h[, .(snp_id, beta_y, se_y, p_y)], by = "snp_id")
    res <- x[, {
      a <- run_coloc(.SD, o, 5e-6)
      b <- run_coloc(.SD, o, 1e-5)
      lx <- which.max(abs(beta / se)); ly <- which.min(p_y)
      c(as.list(a$pp),
        list(PP.H4_p12_1e5 = b$pp[["PP.H4.abf"]], top_snp_H4 = a$top, top_snp_PP = a$top_pp,
             lead_mqtl_snp = snp_id[lx], lead_mqtl_z = beta[lx] / se[lx],
             lead_gwas_snp = snp_id[ly], lead_gwas_p = p_y[ly],
             gwas_p_at_lead_mqtl = p_y[lx]))
    }, by = cpg_id]
    setnames(res, c("nsnps", paste0("PP.H", 0:4, ".abf")), c("nsnps", paste0("PP.H", 0:4)))
    res <- merge(hits[, .(cpg_id, cpg_chr, cpg_pos, n_iv, b, se, p, fdr, or_per_0.1, flag_region)], res,
                 by = "cpg_id", all.x = TRUE)
    res[, coloc := fcase(is.na(PP.H4), "not_tested",
                         PP.H4 >= 0.8, "strong_H4",
                         PP.H4 >= 0.5, "moderate_H4",
                         PP.H3 >= 0.5, "distinct_H3",
                         default = "inconclusive")]
    setorder(res, p)
    res[, `:=`(exposure = ds, outcome = oc)]
    fwrite(res, out_path(cfg, "coloc", sprintf("%s__%s.tsv.gz", ds, oc)), sep = "\t")
    summ[[paste(ds, oc)]] <- data.table(exposure = ds, outcome = oc, mr_hits = nrow(hits),
                                        tested = res[!is.na(PP.H4), .N],
                                        strong_H4 = res[coloc == "strong_H4", .N],
                                        moderate_H4 = res[coloc == "moderate_H4", .N],
                                        distinct_H3 = res[coloc == "distinct_H3", .N],
                                        inconclusive = res[coloc == "inconclusive", .N],
                                        median_nsnps = median(res$nsnps, na.rm = TRUE))
    logmsg("  %s: %d hits, %d tested; H4>=0.8: %d, H3>=0.5: %d", ds, nrow(hits),
           summ[[paste(ds, oc)]]$tested, summ[[paste(ds, oc)]]$strong_H4, summ[[paste(ds, oc)]]$distinct_H3)
  }
  rm(g); gc(verbose = FALSE)
}
s <- rbindlist(summ)
fwrite(s, out_path(cfg, "coloc", "summary.tsv"), sep = "\t")
print(s)
