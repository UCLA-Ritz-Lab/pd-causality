#!/usr/bin/env Rscript
# Phase 1, step 7: coloc-SuSiE for the primary MR hits (allows several causal variants
# per trait, which single-variant coloc.abf cannot).
#
# Why: at multi-signal PD loci (SNCA, TMEM175/GAK, 16p11.2, MAPT ...) coloc.abf compares
# the mQTL with the *strongest* GWAS signal and returns H3 even when the methylation
# signal matches a secondary GWAS signal.
#
# LD set-up (decided 2026-11-15):
#   mQTL  : PEG1 controls regional statistics (step 5) with IN-SAMPLE LD from the same
#           231 individuals (correlation of covariate-residualised ALT dosages, i.e. the
#           LD the regression actually saw). Single cohort rather than the controls meta,
#           so that statistics and LD come from the same sample.
#   GWAS  : Nalls z-scores with 1000G phase 3 EUR LD (n = 503), sign-aligned to ALT.
#   SNPs  : PEG1-controls regional SNPs ∩ harmonised GWAS ∩ EUR panel (biallelic, global
#           MAF >= 5%, polymorphic in EUR).
# SuSiE: coloc::runsusie (L = 10), then coloc.susie (p1 = p2 = 1e-4, p12 = 5e-6). A trait
# with no 95% credible set cannot be tested by coloc.susie; then only the coloc.abf
# result applies. LD mismatch between 1000G and the Nalls meta-analysis can create
# spurious GWAS credible sets; the z-vs-LD consistency check (kriging_rss lambda) is
# reported per region.
#
# Targets: CpGs with FDR < cfg$fdr in results/phase1/mr/meta_controls__nalls2019.tsv.gz.
# Output: results/phase1/coloc_susie/meta_controls__{outcome}.tsv.gz (one row per CpG, best
#         credible-set pair) and *_pairs.tsv.gz (all credible-set pairs).
# Usage:  Rscript scripts/R/07_coloc_susie.R [outcome]   (default nalls2019)

source("scripts/R/phase1_lib.R")
suppressPackageStartupMessages({ library(coloc); library(susieR) })
cfg <- load_config()
args <- commandArgs(trailingOnly = TRUE)
oc <- if (length(args)) args[1] else "nalls2019"
o  <- cfg$outcomes[[oc]]
co <- "peg1_controls"
cc <- cfg$cohorts[[co]]
pipe <- cc$pipeline
set.seed(20261115)

hits <- fread(out_path(cfg, "mr", sprintf("meta_controls__%s.tsv.gz", oc)),
              colClasses = list(character = "cpg_chr"))[fdr < cfg$fdr]
reg <- fread(out_path(cfg, "regional_mqtl", paste0(co, ".tsv.gz")),
             colClasses = list(character = c("cpg_id", "snp_id")))[cpg_id %in% hits$cpg_id]
logmsg("%s: %d hit CpGs; %d with PEG1-controls regional stats", oc, nrow(hits), uniqueN(reg$cpg_id))

# ── alleles + GWAS harmonised to ALT ─────────────────────────────────────────
ea <- fread(pdqtl_path(cfg, pipe, "positions/export_allele.txt"), header = FALSE, col.names = c("id", "counted"))
ea[, snp_id := sub(":[ACGT]+:[ACGT]+$", "", id)]
ea <- ea[snp_id %in% reg$snp_id]
ea[, c("chr", "pos", "ref", "alt") := tstrsplit(id, ":", fixed = TRUE)]
ea[, pos := as.integer(pos)]
reg <- merge(reg, ea[, .(snp_id, chr, pos, ref, alt)], by = "snp_id")
xs <- unique(reg[, .(snp_id, chr, pos, ref, alt, eaf)], by = "snp_id")
g <- read_gwas(o$file)
h <- harmonise_gwas(xs, g[unique(xs[, .(chr, pos)]), on = c("chr", "pos"), nomatch = NULL],
                    cfg$instruments)[is.na(drop)]
rm(g); gc(verbose = FALSE)
reg <- merge(reg, h[, .(snp_id, beta_y, se_y)], by = "snp_id")

# ── PEG1 controls: covariate-residualised dosages (individual level, in memory only) ──
ddir <- pdqtl_path(cfg, pipe, "datasets", co)
read_header <- function(f, sep) { con <- gzcon(file(f, "rb")); on.exit(close(con)); trimws(strsplit(readLines(con, n = 1), sep)[[1]]) }
geno_f <- pdqtl_path(cfg, pipe, "genotypes", paste0(co, ".geno.tsv.gz"))
cv <- read.csv(file.path(ddir, "covariates.csv"), row.names = 1, check.names = FALSE)
colnames(cv) <- sub("^X", "", trimws(colnames(cv)))
idm <- read.csv(file.path(ddir, "id_map.csv"), colClasses = "character")
mh <- read_header(file.path(ddir, "methylation.csv.gz"), ",")[-1]
gh <- read_header(geno_f, "\t")[-1]
al <- idm[idm[[1]] %in% mh & idm[[1]] %in% colnames(cv) & idm[[2]] %in% gh, ]
stopifnot(nrow(al) == cc$n_samples)
C <- cbind(1, t(as.matrix(cv[, al[[1]], drop = FALSE]))); storage.mode(C) <- "double"
Qc <- qr.Q(qr(C))
kf <- tempfile(); writeLines(unique(reg$snp_id), kf)
gd <- fread(cmd = sprintf("zcat %s | awk -F'\\t' 'NR==FNR{k[$1]=1; next} FNR==1 || ($1 in k)' %s -", shQuote(geno_f), kf),
            header = TRUE, sep = "\t", colClasses = list(character = 1))
setnames(gd, read_header(geno_f, "\t"))
stopifnot(!anyDuplicated(gd[[1]]))
Gp <- as.matrix(gd[, al[[2]], with = FALSE]); rownames(Gp) <- gd[[1]]; rm(gd)
Gp <- t(Gp)
Gp <- Gp - Qc %*% crossprod(Qc, Gp)
logmsg("PEG1 controls: %d samples, %s SNPs loaded", nrow(Gp), format(ncol(Gp), big.mark = ","))

# ── 1000G EUR dosages, aligned to ALT ───────────────────────────────────────
panel <- fread(cfg$clump$panel)
keep_f <- tempfile(); writeLines(c("#IID", panel[super_pop == cfg$clump$super_pop, sample]), keep_f)
eur <- list()
for (ch in unique(reg$chr)) {
  pref <- sub("\\{chrom\\}", ch, cfg$clump$ref_prefix)
  pv <- fread(cmd = sprintf("grep -v '^##' %s.pvar | cut -f1-5", pref), colClasses = list(character = c("#CHROM", "ID")))
  setnames(pv, c("chr", "pos", "kid", "kref", "kalt"))
  s <- unique(reg[chr == ch, .(snp_id, pos, ref, alt)])
  m <- merge(s, pv[, !"chr"], by = "pos")[(ref == kref & alt == kalt) | (ref == kalt & alt == kref)][!duplicated(snp_id)]
  idf <- tempfile(); writeLines(m$kid, idf); outp <- tempfile()
  stopifnot(system(sprintf("%s --pfile %s --keep %s --extract %s --export Av --out %s --threads 4 --memory 8000 > /dev/null 2>&1",
                           cfg$clump$plink2, pref, keep_f, idf, outp)) == 0)
  tr <- fread(paste0(outp, ".traw"))
  X <- as.matrix(tr[, -(1:6)])
  ids <- m$snp_id[match(tr$SNP, m$kid)]
  alt <- m$alt[match(tr$SNP, m$kid)]
  X[tr$COUNTED != alt, ] <- 2 - X[tr$COUNTED != alt, ]          # count ALT
  rownames(X) <- ids
  eur[[ch]] <- X
}
logmsg("1000G EUR dosages loaded for %s SNPs", format(sum(vapply(eur, nrow, 0L)), big.mark = ","))

ld_from <- function(M) { R <- suppressWarnings(cor(M)); R[is.na(R)] <- 0; dimnames(R) <- list(colnames(M), colnames(M)); R }
gp_poly <- colnames(Gp)[apply(Gp, 2, sd) > 1e-8]

# ── per CpG ──────────────────────────────────────────────────────────────────
res <- list(); pairs_out <- list()
only <- Sys.getenv("SUSIE_ONLY")
cpgs <- unique(reg$cpg_id)
if (nzchar(only)) cpgs <- intersect(cpgs, strsplit(only, ",")[[1]])
for (cp in cpgs) {
  d <- reg[cpg_id == cp]
  ch <- d$chr[1]
  Xe <- eur[[ch]]
  sd_e <- apply(Xe[intersect(d$snp_id, rownames(Xe)), , drop = FALSE], 1, sd)
  keep <- intersect(intersect(d$snp_id, names(sd_e)[sd_e > 0]), gp_poly)
  d <- d[snp_id %in% keep]
  setorder(d, pos)
  R1 <- ld_from(Gp[, d$snp_id, drop = FALSE])
  R2 <- ld_from(t(Xe[d$snp_id, , drop = FALSE]))
  maf <- pmin(d$eaf, 1 - d$eaf)
  ds1 <- list(beta = setNames(d$beta, d$snp_id), varbeta = setNames(d$se^2, d$snp_id), snp = d$snp_id,
              type = "quant", N = cc$n_samples, MAF = setNames(maf, d$snp_id), LD = R1)
  ds2 <- list(beta = setNames(d$beta_y, d$snp_id), varbeta = setNames(d$se_y^2, d$snp_id), snp = d$snp_id,
              type = "cc", s = o$n_cases / (o$n_cases + o$n_controls), N = o$n_cases + o$n_controls,
              MAF = setNames(maf, d$snp_id), LD = R2)
  z1 <- d$beta / d$se; z2 <- d$beta_y / d$se_y
  lam1 <- tryCatch(estimate_s_rss(z1, R1, n = cc$n_samples), error = function(e) NA_real_)
  lam2 <- tryCatch(estimate_s_rss(z2, R2, n = ds2$N), error = function(e) NA_real_)
  s1 <- tryCatch(suppressWarnings(suppressMessages(runsusie(ds1, L = 10, repeat_until_convergence = FALSE))), error = function(e) NULL)
  s2 <- tryCatch(suppressWarnings(suppressMessages(runsusie(ds2, L = 10, repeat_until_convergence = FALSE))), error = function(e) NULL)
  ncs1 <- if (is.null(s1)) NA_integer_ else length(s1$sets$cs)
  ncs2 <- if (is.null(s2)) NA_integer_ else length(s2$sets$cs)
  row <- data.table(cpg_id = cp, nsnps = nrow(d), cs_mqtl = ncs1, cs_gwas = ncs2,
                    s_lambda_mqtl = lam1, s_lambda_gwas = lam2,
                    best_H4 = NA_real_, best_H3 = NA_real_, hit_mqtl = NA_character_, hit_gwas = NA_character_)
  if (isTRUE(ncs1 > 0) && isTRUE(ncs2 > 0)) {
    cs <- suppressWarnings(suppressMessages(coloc.susie(s1, s2, p1 = 1e-4, p2 = 1e-4, p12 = 5e-6)))
    sm <- as.data.table(cs$summary)
    if (nrow(sm)) {
      sm[, cpg_id := cp]
      pairs_out[[cp]] <- sm
      b <- sm[which.max(PP.H4.abf)]
      row[, `:=`(best_H4 = b$PP.H4.abf, best_H3 = b$PP.H3.abf, hit_mqtl = b$hit1, hit_gwas = b$hit2)]
    }
  }
  res[[cp]] <- row
  logmsg("  %s: %d SNPs, CS mQTL %s, CS GWAS %s, best H4 %s", cp, nrow(d), ncs1, ncs2,
         ifelse(is.na(row$best_H4), "NA", sprintf("%.2f", row$best_H4)))
}

r <- rbindlist(res)
abf <- fread(out_path(cfg, "coloc", sprintf("meta_controls__%s.tsv.gz", oc)))[, .(cpg_id, abf_H3 = PP.H3, abf_H4 = PP.H4, abf_call = coloc)]
r <- merge(hits[, .(cpg_id, cpg_chr, cpg_pos, b, p, fdr, flag_region)], r, by = "cpg_id", all.x = TRUE)
r <- merge(r, abf, by = "cpg_id", all.x = TRUE)
r[, susie_call := fcase(is.na(cs_mqtl) | is.na(cs_gwas), "not_run",
                        cs_mqtl == 0 | cs_gwas == 0, "no_credible_set",
                        best_H4 >= 0.8, "strong_H4",
                        best_H4 >= 0.5, "moderate_H4",
                        best_H3 >= 0.5, "distinct_H3",
                        default = "inconclusive")]
setorder(r, p)
fwrite(r, out_path(cfg, "coloc_susie", sprintf("meta_controls__%s.tsv.gz", oc)), sep = "\t")
if (length(pairs_out)) fwrite(rbindlist(pairs_out, fill = TRUE),
                              out_path(cfg, "coloc_susie", sprintf("meta_controls__%s_pairs.tsv.gz", oc)), sep = "\t")
print(r[, .N, by = .(abf_call, susie_call)][order(abf_call, susie_call)])
