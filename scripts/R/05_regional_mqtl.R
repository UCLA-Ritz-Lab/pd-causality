#!/usr/bin/env Rscript
# Phase 1, step 5: untruncated regional cis-mQTL statistics for the MR hit CpGs (for coloc).
#
# The PD-QTLs cis files keep only p < 1e-3 (PPMI) / 1e-2 (PEG), which biases coloc, so
# this refits every SNP within +-1 Mb of each target CpG using exactly the PD-QTLs model:
# methylation beta-value ~ ALT dosage + covariates (MatrixEQTL modelLINEAR), on the same
# prepared inputs (PD-QTLs results/{PIPE}/datasets/{cohort}/ methylation.csv.gz,
# covariates.csv, id_map.csv; genotypes/{cohort}.geno.tsv.gz; SNP universe
# positions/snp_positions.csv). Fit by residualising methylation and dosage on
# [intercept + covariates] (Frisch-Waugh), as MatrixEQTL does:
#   beta = <g~, y~> / <g~, g~>,  r = cor(g~, y~),  df = n - 2 - n_cov,
#   t = r sqrt(df / (1 - r^2)),  se = |beta / t|.
# Validated against the existing truncated cis file on shared pairs (must agree).
#
# Targets: CpGs with FDR < cfg$fdr in any results/phase1/mr/*__nalls2019.tsv.gz.
# Output: results/phase1/regional_mqtl/{cohort}.tsv.gz (summary statistics only) and
#         validation.tsv. No individual-level values are written or printed.
# Usage:  Rscript scripts/R/05_regional_mqtl.R [cohort ...]   (default: all in config)

source("scripts/R/phase1_lib.R")
cfg  <- load_config()
args <- commandArgs(trailingOnly = TRUE)
cohorts <- if (length(args)) args else names(cfg$cohorts)
WIN <- 1e6

mr_files <- list.files(out_path(cfg, "mr"), pattern = "__nalls2019\\.tsv\\.gz$", full.names = TRUE)
targets <- unique(rbindlist(lapply(mr_files, function(f) fread(f, select = c("cpg_id", "fdr"))[fdr < cfg$fdr]))$cpg_id)
logmsg("target CpGs (FDR < %.2f in any nalls2019 analysis): %d", cfg$fdr, length(targets))

read_header <- function(f, sep) {
  con <- gzcon(file(f, "rb")); on.exit(close(con))
  trimws(strsplit(readLines(con, n = 1), sep)[[1]])
}
# rows of a gzipped table whose first field is in `keys` (plus the header)
read_rows <- function(f, keys, sep) {
  kf <- tempfile(); writeLines(keys, kf); on.exit(unlink(kf))
  d <- fread(cmd = sprintf("zcat %s | awk -F'%s' 'NR==FNR{k[$1]=1; next} FNR==1 || ($1 in k)' %s -",
                           shQuote(f), sep, kf),
             header = TRUE, sep = sep, colClasses = list(character = 1))
  setnames(d, read_header(f, sep))
  d
}

val_rows <- list()
for (co in cohorts) {
  cc <- cfg$cohorts[[co]]
  pipe <- cc$pipeline
  ddir <- pdqtl_path(cfg, pipe, "datasets", co)
  logmsg("=== %s (%s)", co, pipe)

  # ── samples: methylation x covariates x genotypes ──────────────────────────
  meth_f <- file.path(ddir, "methylation.csv.gz")
  geno_f <- pdqtl_path(cfg, pipe, "genotypes", paste0(co, ".geno.tsv.gz"))
  cv <- read.csv(file.path(ddir, "covariates.csv"), row.names = 1, check.names = FALSE)
  colnames(cv) <- sub("^X", "", trimws(colnames(cv)))
  idm <- read.csv(file.path(ddir, "id_map.csv"), colClasses = "character")
  meth_id <- names(idm)[1]; geno_id <- names(idm)[2]   # PEG: Pegid,GWAS_ID ; PPMI: PATNO,plink_iid
  mh <- read_header(meth_f, ",")[-1]
  gh <- read_header(geno_f, "\t")[-1]
  al <- idm[idm[[meth_id]] %in% mh & idm[[meth_id]] %in% colnames(cv) & idm[[geno_id]] %in% gh, ]
  n <- nrow(al); n_cov <- nrow(cv)
  stopifnot(n_cov == cc$n_cov, !anyDuplicated(al[[meth_id]]), !anyDuplicated(al[[geno_id]]))
  df <- n - 2 - n_cov
  logmsg("  aligned samples: %d (config n_samples %d); covariates: %d", n, cc$n_samples, n_cov)

  C <- cbind(1, t(as.matrix(cv[, al[[meth_id]], drop = FALSE])))
  storage.mode(C) <- "double"
  Qc <- qr.Q(qr(C))
  resid <- function(M) M - Qc %*% crossprod(Qc, M)

  # ── methylation for targets ────────────────────────────────────────────────
  m <- read_rows(meth_f, targets, ",")
  Y <- t(as.matrix(m[, al[[meth_id]], with = FALSE])); colnames(Y) <- m[[1]]
  stopifnot(!anyNA(Y))
  Yr <- resid(Y)
  logmsg("  target CpGs on this array: %d of %d", ncol(Y), length(targets))

  # ── window SNPs ────────────────────────────────────────────────────────────
  cpos <- fread(pdqtl_path(cfg, cc$probes), colClasses = "character")[cpg_id %in% colnames(Y)]
  cpos[, `:=`(lo = as.integer(pos) - WIN, hi = as.integer(pos) + WIN)]
  spos <- fread(pdqtl_path(cfg, pipe, "positions", "snp_positions.csv"), colClasses = "character")
  setnames(spos, 1:3, c("snp_id", "chr", "pos"))
  spos[, pos := as.integer(pos)]
  # PPMI's snp_positions.csv repeats ~2k chr:pos ids (multiallelic sites); the genotype
  # matrix has one record per position, so one row per id is all that can be tested.
  spos <- unique(spos, by = "snp_id")
  pairs <- spos[cpos, on = .(chr, pos >= lo, pos <= hi), nomatch = NULL,
                .(cpg_id = i.cpg_id, snp_id = x.snp_id, snp_pos = x.pos, cpg_pos = as.integer(i.pos))]
  snps <- unique(pairs$snp_id)
  logmsg("  window SNPs: %s; SNP-CpG pairs: %s", format(length(snps), big.mark = ","),
         format(nrow(pairs), big.mark = ","))
  g <- read_rows(geno_f, snps, "\t")
  stopifnot(!anyDuplicated(g[[1]]))
  G <- as.matrix(g[, al[[geno_id]], with = FALSE]); rownames(G) <- g[[1]]
  rm(g); gc(verbose = FALSE)
  G <- t(G)                                                    # samples x SNPs
  if (anyNA(G)) { mu <- colMeans(G, na.rm = TRUE); idx <- which(is.na(G), arr.ind = TRUE); G[idx] <- mu[idx[, 2]] }
  alt_freq <- colMeans(G) / 2
  Gr <- resid(G)
  gss <- colSums(Gr^2)

  # ── per-pair stats, CpG by CpG ─────────────────────────────────────────────
  out <- pairs[snp_id %in% colnames(G)][, {
    gi <- match(snp_id, colnames(G))
    y <- Yr[, cpg_id[1]]
    gy <- as.vector(crossprod(Gr[, gi, drop = FALSE], y))
    b <- gy / gss[gi]
    r <- gy / sqrt(gss[gi] * sum(y^2))
    tt <- r * sqrt(df / (1 - r^2))
    .(snp_id, dist_bp = snp_pos - cpg_pos, beta = b, se = abs(b / tt), t = tt,
      p = 2 * pt(-abs(tt), df), eaf = alt_freq[gi])
  }, by = cpg_id]
  out <- out[is.finite(t)]
  out[, `:=`(cohort = co, n = n)]
  fwrite(out, out_path(cfg, "regional_mqtl", paste0(co, ".tsv.gz")), sep = "\t")
  logmsg("  wrote %s pair statistics", format(nrow(out), big.mark = ","))

  # ── validate against PD-QTLs cis file (shared pairs) ───────────────────────
  cis_f <- pdqtl_path(cfg, pipe, "mqtl", co, "cis_mqtls.txt.gz")
  cf <- tempfile(); writeLines(colnames(Y), cf)
  cis <- fread(cmd = sprintf("zcat %s | awk -F'\\t' 'NR==FNR{k[$1]=1; next} FNR==1 || ($2 in k)' %s -",
                             shQuote(cis_f), cf), colClasses = list(character = 1:2))
  cis[, SNP := sub(":[ACGT]+:[ACGT]+$", "", SNP)]
  v <- merge(out, cis, by.x = c("cpg_id", "snp_id"), by.y = c("gene", "SNP"), suffixes = c("", ".cis"))
  val_rows[[co]] <- data.table(cohort = co, n = n, shared_pairs = nrow(v),
                               cor_beta = cor(v$beta, v$beta.cis), max_abs_diff_beta = max(abs(v$beta - v$beta.cis)),
                               cor_t = cor(v$t, v$`t-stat`), max_rel_diff_t = max(abs(v$t / v$`t-stat` - 1)))
  print(val_rows[[co]])
  stopifnot(nrow(v) > 0, val_rows[[co]]$max_rel_diff_t < 1e-3)
}
vf <- out_path(cfg, "regional_mqtl", "validation.tsv")
vv <- rbindlist(val_rows)
if (file.exists(vf) && length(args)) vv <- rbind(fread(vf)[!cohort %in% vv$cohort], vv)
fwrite(vv, vf, sep = "\t")
