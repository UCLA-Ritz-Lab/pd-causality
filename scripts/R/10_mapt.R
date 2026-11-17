#!/usr/bin/env Rscript
# Phase 1, step 10: MAPT 17q21.31 H1/H2 handling.
#
# Problem. The ~900 kb H1/H2 inversion keeps a ~1.5 Mb block in near-complete LD. H2 is
# protective for PD and changes methylation (and expression) of many CpGs/genes in the
# block. A single-SNP Wald ratio for a CpG whose instrument tags H2 therefore estimates
# "effect of the H2 haplotype on PD, scaled by its effect on this CpG". Every H2-regulated
# CpG gets a significant estimate, so it cannot point to a CpG (exclusion restriction
# fails by construction: one instrument, many mediators).
#
# What this step does, for MAPT-region CpGs with MR FDR < cfg$fdr in any *__nalls2019
# analysis (those with regional stats from step 5):
#   1. H2 tag: rs8070723, 17:44081064 A>G, G (ALT) = H2. Reports how strongly each CpG's
#      primary instrument tags H2 (r^2, 1000G EUR) and the haplotype-level mQTL/Wald ratio.
#   2. Conditional mQTL (exact, in-sample): refits the regional mQTL in PEG1 and PPMI
#      controls with H2-tag dosage added to the covariates (model otherwise as step 5),
#      then common-effect IVW meta over SNPs with identical alleles in both cohorts.
#   3. Conditional GWAS (approximate, summary-level): for SNP j and tag k with signed LD
#      r (1000G EUR, ALT-coded), z_j|k = (z_j - r z_k) / sqrt(1 - r^2),
#      se_j|k = se_j / sqrt(1 - r^2). Assumes 1000G EUR LD matches the GWAS sample.
#   4. H2-independent signal: if a conditional mQTL at p < 5e-8 remains, conditional MR
#      (Wald ratio at the conditional mQTL lead SNP) and conditional coloc.abf
#      (p12 = 5e-6) on the conditional statistics.
#   5. H2 effect in every PD GWAS (Nalls, no-proxy, in-house) and the number of CpGs
#      genome-wide whose cis-mQTL includes the H2 tag at p < 5e-8 (breadth of pleiotropy).
# Classification per CpG:
#   H2_haplotype       : instrument (or, for non-primary hits, controls-meta lead SNP)
#                        r^2 with H2 >= 0.5, and no conditional mQTL at p < 5e-8
#   H2_tagging_residual: instrument tags H2 (r^2 >= 0.5) but a conditional mQTL remains
#   not_H2             : instrument r^2 with H2 < 0.5
# Summary statistics only. Output: results/phase1/mapt/{cpgs,gwas_h2,breadth}.tsv
# Usage:  Rscript scripts/R/10_mapt.R

source("scripts/R/phase1_lib.R")
suppressPackageStartupMessages(library(coloc))
cfg <- load_config()
WIN <- 1e6
H2  <- list(snp_id = "17:44081064", ref = "A", alt = "G", rsid = "rs8070723")
COHORTS <- c("peg1_controls", "ppmi_controls")
P_GW <- 5e-8

# ── Targets ──────────────────────────────────────────────────────────────────
mr_files <- list.files(out_path(cfg, "mr"), pattern = "__nalls2019\\.tsv\\.gz$", full.names = TRUE)
mr_all <- rbindlist(lapply(mr_files, function(f) fread(f, select = c("cpg_id", "cpg_chr", "cpg_pos", "fdr", "flag_region"),
                                                      colClasses = list(character = "cpg_chr"))[, analysis := sub("\\.tsv\\.gz$", "", basename(f))]))
targets <- unique(mr_all[fdr < cfg$fdr & flag_region == "MAPT_17q21.31", .(cpg_id, cpg_pos)])
logmsg("MAPT-region CpGs with FDR < %.2f in any nalls2019 analysis: %d", cfg$fdr, nrow(targets))

read_header <- function(f, sep) { con <- gzcon(file(f, "rb")); on.exit(close(con)); trimws(strsplit(readLines(con, n = 1), sep)[[1]]) }
read_rows <- function(f, keys, sep) {
  kf <- tempfile(); writeLines(keys, kf); on.exit(unlink(kf))
  d <- fread(cmd = sprintf("zcat %s | awk -F'%s' 'NR==FNR{k[$1]=1; next} FNR==1 || ($1 in k)' %s -", shQuote(f), sep, kf),
             header = TRUE, sep = sep, colClasses = list(character = 1))
  setnames(d, read_header(f, sep)); d
}

# ── 2. Marginal + H2-conditional regional mQTL per control cohort ────────────
fit_pairs <- function(pairs, Gr, Yr, df) {
  gss <- colSums(Gr^2)
  pairs[snp_id %in% colnames(Gr)][, {
    gi <- match(snp_id, colnames(Gr)); y <- Yr[, cpg_id[1]]
    gy <- as.vector(crossprod(Gr[, gi, drop = FALSE], y))
    b <- gy / gss[gi]; r <- gy / sqrt(gss[gi] * sum(y^2)); tt <- r * sqrt(df / (1 - r^2))
    .(snp_id, beta = b, se = abs(b / tt), t = tt)
  }, by = cpg_id][is.finite(t)]
}

reg <- list(); h2freq <- list()
for (co in COHORTS) {
  cc <- cfg$cohorts[[co]]; pipe <- cc$pipeline
  ddir <- pdqtl_path(cfg, pipe, "datasets", co)
  meth_f <- file.path(ddir, "methylation.csv.gz")
  geno_f <- pdqtl_path(cfg, pipe, "genotypes", paste0(co, ".geno.tsv.gz"))
  cv <- read.csv(file.path(ddir, "covariates.csv"), row.names = 1, check.names = FALSE)
  colnames(cv) <- sub("^X", "", trimws(colnames(cv)))
  idm <- read.csv(file.path(ddir, "id_map.csv"), colClasses = "character")
  mid <- names(idm)[1]; gid <- names(idm)[2]
  al <- idm[idm[[mid]] %in% read_header(meth_f, ",")[-1] & idm[[mid]] %in% colnames(cv) &
            idm[[gid]] %in% read_header(geno_f, "\t")[-1], ]
  n <- nrow(al); n_cov <- nrow(cv)
  stopifnot(n == cc$n_samples, n_cov == cc$n_cov)

  m <- read_rows(meth_f, targets$cpg_id, ",")
  Y <- t(as.matrix(m[, al[[mid]], with = FALSE])); colnames(Y) <- m[[1]]
  cpos <- targets[cpg_id %in% colnames(Y)]
  spos <- fread(pdqtl_path(cfg, pipe, "positions", "snp_positions.csv"), colClasses = "character")
  setnames(spos, 1:3, c("snp_id", "chr", "pos"))
  spos <- unique(spos[chr == "17"], by = "snp_id")[, pos := as.integer(pos)]
  cpos[, `:=`(lo = cpg_pos - WIN, hi = cpg_pos + WIN)]
  pairs <- spos[cpos, on = .(pos >= lo, pos <= hi), nomatch = NULL, .(cpg_id = i.cpg_id, snp_id = x.snp_id)]
  g <- read_rows(geno_f, unique(c(pairs$snp_id, H2$snp_id)), "\t")
  stopifnot(!anyDuplicated(g[[1]]), H2$snp_id %in% g[[1]])
  G <- t(as.matrix(g[, al[[gid]], with = FALSE])); colnames(G) <- g[[1]]; rm(g)
  if (anyNA(G)) { mu <- colMeans(G, na.rm = TRUE); ix <- which(is.na(G), arr.ind = TRUE); G[ix] <- mu[ix[, 2]] }
  h2 <- G[, H2$snp_id]
  h2freq[[co]] <- data.table(cohort = co, n = n, h2_freq = mean(h2) / 2)

  C0 <- cbind(1, t(as.matrix(cv[, al[[mid]], drop = FALSE]))); storage.mode(C0) <- "double"
  fits <- list()
  for (cond in c(FALSE, TRUE)) {
    Cm <- if (cond) cbind(C0, h2) else C0
    Q <- qr.Q(qr(Cm)); res <- function(M) M - Q %*% crossprod(Q, M)
    df <- n - 2 - (ncol(Cm) - 1)
    f <- fit_pairs(pairs, res(G), res(Y), df)
    f[, p := 2 * pt(-abs(t), df)]
    if (cond) f <- f[snp_id != H2$snp_id]
    fits[[if (cond) "cond" else "marg"]] <- f
  }
  ea <- fread(pdqtl_path(cfg, pipe, "positions/export_allele.txt"), header = FALSE, col.names = c("id", "counted"))
  ea[, snp_id := sub(":[ACGT]+:[ACGT]+$", "", id)]
  ea <- ea[snp_id %in% colnames(G)][, c("chr", "pos", "ref", "alt") := tstrsplit(id, ":", fixed = TRUE)]
  stopifnot(ea[snp_id == H2$snp_id, ref == H2$ref & alt == H2$alt])
  af <- data.table(snp_id = colnames(G), eaf = colMeans(G) / 2)
  for (k in names(fits)) {
    fits[[k]] <- merge(merge(fits[[k]], ea[, .(snp_id, chr, pos = as.integer(pos), ref, alt)], by = "snp_id"), af, by = "snp_id")
    fits[[k]][, `:=`(cohort = co, n = n)]
  }
  reg[[co]] <- fits
  logmsg("%s: n = %d, %d CpGs, H2 frequency %.3f", co, n, ncol(Y), mean(h2) / 2)
}
# sanity: marginal refit must match step 5's regional file (same model)
for (co in COHORTS) {
  s5 <- fread(out_path(cfg, "regional_mqtl", paste0(co, ".tsv.gz")), colClasses = list(character = c("cpg_id", "snp_id")))
  v <- merge(reg[[co]]$marg, s5[, .(cpg_id, snp_id, t5 = t)], by = c("cpg_id", "snp_id"))
  stopifnot(nrow(v) > 0, max(abs(v$t / v$t5 - 1)) < 1e-6)
}

meta <- function(k) {
  d <- rbindlist(lapply(COHORTS, function(co) reg[[co]][[k]]))
  d[, ok := .N == length(COHORTS) & uniqueN(paste(ref, alt)) == 1, by = .(cpg_id, snp_id)]
  d <- d[ok == TRUE][, w := 1 / se^2]
  d[, .(chr = chr[1], pos = pos[1], ref = ref[1], alt = alt[1], beta = sum(w * beta) / sum(w), se = sqrt(1 / sum(w)),
        eaf = sum(n * eaf) / sum(n), n = sum(n)), by = .(cpg_id, snp_id)][, z := beta / se]
}
mx <- meta("marg"); cx <- meta("cond")

# ── 3. GWAS (Nalls) harmonised to ALT + EUR LD with the H2 tag ────────────────
xs <- unique(mx[, .(snp_id, chr, pos, ref, alt, eaf)])
g <- read_gwas(cfg$outcomes$nalls2019$file)
gy <- harmonise_gwas(xs, g[unique(xs[, .(chr, pos)]), on = c("chr", "pos"), nomatch = NULL], cfg$instruments)[is.na(drop)]
rm(g); gc(verbose = FALSE)
gy <- gy[, .(snp_id, beta_y, se_y, z_y = beta_y / se_y, p_y)]
stopifnot(H2$snp_id %in% gy$snp_id)

panel <- fread(cfg$clump$panel)
keep_f <- tempfile(); writeLines(c("#IID", panel[super_pop == cfg$clump$super_pop, sample]), keep_f)
pref <- sub("\\{chrom\\}", "17", cfg$clump$ref_prefix)
pv <- fread(cmd = sprintf("grep -v '^##' %s.pvar | cut -f1-5", pref), colClasses = list(character = c("#CHROM", "ID")))
setnames(pv, c("chr", "pos", "kid", "kref", "kalt"))
mm <- merge(xs[, .(snp_id, pos, ref, alt)], pv[, !"chr"], by = "pos")[(ref == kref & alt == kalt) | (ref == kalt & alt == kref)][!duplicated(snp_id)]
idf <- tempfile(); writeLines(mm$kid, idf); outp <- tempfile()
stopifnot(system(sprintf("%s --pfile %s --keep %s --extract %s --export Av --out %s --threads 4 --memory 8000 > /dev/null 2>&1",
                         cfg$clump$plink2, pref, keep_f, idf, outp)) == 0)
tr <- fread(paste0(outp, ".traw"))
X <- as.matrix(tr[, -(1:6)])
alt <- mm$alt[match(tr$SNP, mm$kid)]
X[tr$COUNTED != alt, ] <- 2 - X[tr$COUNTED != alt, ]
rownames(X) <- mm$snp_id[match(tr$SNP, mm$kid)]
X <- X[apply(X, 1, sd) > 0, , drop = FALSE]
stopifnot(H2$snp_id %in% rownames(X))
ld_h2 <- data.table(snp_id = rownames(X), r_h2 = as.vector(cor(t(X), X[H2$snp_id, ])))
logmsg("EUR LD with H2 tag: %s SNPs; %d with r^2 >= 0.8", format(nrow(ld_h2), big.mark = ","), ld_h2[r_h2^2 >= 0.8, .N])

zk <- gy[snp_id == H2$snp_id, z_y]
gc_ <- merge(gy, ld_h2, by = "snp_id")[snp_id != H2$snp_id & r_h2^2 < 0.95]
gc_[, `:=`(z_y_c = (z_y - r_h2 * zk) / sqrt(1 - r_h2^2), se_y_c = se_y / sqrt(1 - r_h2^2))]
gc_[, `:=`(beta_y_c = z_y_c * se_y_c, p_y_c = 2 * pnorm(-abs(z_y_c)))]

# ── 1 + 4. Per-CpG summary ───────────────────────────────────────────────────
inst <- fread(out_path(cfg, "instruments", "meta_controls__nalls2019.tsv.gz"), colClasses = list(character = c("snp_id", "cpg_id")))
inst <- inst[cpg_id %in% targets$cpg_id][order(-neglog10p)][!duplicated(cpg_id), .(cpg_id, iv_snp = snp_id)]
prim <- fread(out_path(cfg, "mr", "meta_controls__nalls2019.tsv.gz"))[cpg_id %in% targets$cpg_id, .(cpg_id, primary_b = b, primary_p = p, primary_fdr = fdr)]
h2y <- gy[snp_id == H2$snp_id]

run_abf <- function(d) {
  o <- cfg$outcomes$nalls2019; maf <- pmin(d$eaf, 1 - d$eaf)
  capture.output(r <- suppressWarnings(suppressMessages(coloc.abf(
    list(beta = d$beta, varbeta = d$se^2, snp = d$snp_id, type = "quant", N = d$n[1], MAF = maf),
    list(beta = d$beta_y_c, varbeta = d$se_y_c^2, snp = d$snp_id, type = "cc",
         s = o$n_cases / (o$n_cases + o$n_controls), N = o$n_cases + o$n_controls, MAF = maf),
    p1 = 1e-4, p2 = 1e-4, p12 = 5e-6))))
  r$summary
}

rows <- list()
for (cp in targets$cpg_id) {
  if (!cp %in% mx$cpg_id) next
  m0 <- mx[cpg_id == cp]; c0 <- cx[cpg_id == cp]
  mh <- m0[snp_id == H2$snp_id]
  iv <- inst[cpg_id == cp, iv_snp]
  iv_r2 <- if (length(iv) && iv %in% ld_h2$snp_id) ld_h2[snp_id == iv, r_h2^2] else NA_real_
  ml <- m0[which.max(abs(z))]
  cl <- c0[which.max(abs(z))]
  # conditional MR at the conditional mQTL lead SNP present in the conditional GWAS
  cg <- merge(c0, gc_[, .(snp_id, beta_y_c, se_y_c, p_y_c)], by = "snp_id")
  cgl <- cg[which.max(abs(z))]
  pp <- if (nrow(cg) >= 50 && min(2 * pnorm(-abs(cg$z))) < P_GW) run_abf(cg) else NULL
  rows[[cp]] <- data.table(
    cpg_id = cp, cpg_pos = targets[cpg_id == cp, cpg_pos], iv_snp = if (length(iv)) iv else NA_character_,
    iv_r2_h2_eur = iv_r2,
    marg_lead_snp = ml$snp_id, marg_lead_z = ml$z, marg_lead_r2_h2 = ld_h2[snp_id == ml$snp_id, r_h2^2][1],
    h2_mqtl_beta = mh$beta, h2_mqtl_z = mh$z,
    h2_wald_b = h2y$beta_y / mh$beta, h2_wald_se = h2y$se_y / abs(mh$beta),
    cond_lead_snp = cl$snp_id, cond_lead_z = cl$z, cond_lead_p = 2 * pnorm(-abs(cl$z)),
    cond_n_gw = c0[2 * pnorm(-abs(z)) < P_GW, .N],
    cond_mr_snp = cgl$snp_id, cond_mr_bx = cgl$beta, cond_mr_px = 2 * pnorm(-abs(cgl$z)),
    cond_mr_b = cgl$beta_y_c / cgl$beta, cond_mr_se = cgl$se_y_c / abs(cgl$beta), cond_gwas_p = cgl$p_y_c,
    cond_abf_H3 = if (is.null(pp)) NA_real_ else pp[["PP.H3.abf"]],
    cond_abf_H4 = if (is.null(pp)) NA_real_ else pp[["PP.H4.abf"]])
}
cpgs <- rbindlist(rows)
cpgs[, cond_mr_p := 2 * pnorm(-abs(cond_mr_b / cond_mr_se))]
cpgs[cond_mr_px >= P_GW, c("cond_mr_b", "cond_mr_se", "cond_mr_p") := NA_real_]
# CpGs that are hits only in non-primary analyses have no primary instrument: use the
# marginal controls-meta lead SNP instead.
cpgs[, r2_h2 := fifelse(is.na(iv_r2_h2_eur), marg_lead_r2_h2, iv_r2_h2_eur)]
cpgs[, mapt_class := fcase(is.na(r2_h2), "lead_not_in_ref",
                           r2_h2 >= 0.5 & cond_n_gw == 0, "H2_haplotype",
                           r2_h2 >= 0.5, "H2_tagging_residual",
                           default = "not_H2")]
cpgs <- merge(cpgs, prim, by = "cpg_id", all.x = TRUE)
setorder(cpgs, cpg_pos)
fwrite(cpgs, out_path(cfg, "mapt", "cpgs.tsv"), sep = "\t")

# ── 5a. H2 effect in each PD GWAS ────────────────────────────────────────────
gw <- list()
for (oc in names(cfg$outcomes)) {
  o <- cfg$outcomes[[oc]]
  r <- read_gwas(o$file)[chr == "17" & pos == 44081064L]
  r <- r[(ea == H2$alt & oa == H2$ref) | (ea == H2$ref & oa == H2$alt)]
  if (!nrow(r)) next
  s <- fifelse(r$ea == H2$alt, 1, -1)
  gw[[oc]] <- data.table(outcome = oc, n_cases = o$n_cases, n_controls = o$n_controls,
                         h2_logOR = s * r$beta_g, se = r$se_y, p = r$p_y,
                         h2_freq = fifelse(s > 0, r$eaf_g, 1 - r$eaf_g))
}
gw <- rbindlist(gw)
fwrite(gw, out_path(cfg, "mapt", "gwas_h2.tsv"), sep = "\t")

# ── 5b. Breadth: CpGs whose cis-mQTL includes the H2 tag at p < 5e-8 ─────────
br <- list()
for (co in c(COHORTS, "peg1_cases", "ppmi_cases")) {
  pipe <- cfg$cohorts[[co]]$pipeline
  f <- pdqtl_path(cfg, pipe, "mqtl", co, "cis_mqtls.txt.gz")
  d <- fread(cmd = sprintf("zcat %s | awk -F'\\t' 'NR==1 || $1==\"%s\" || index($1, \"%s:\")==1'", shQuote(f), H2$snp_id, H2$snp_id))
  pc <- grep("^p", names(d), value = TRUE)[1]
  br[[co]] <- data.table(cohort = co, cpgs_h2_p5e8 = d[get(pc) < P_GW, uniqueN(gene)],
                         cpgs_h2_p1e5 = d[get(pc) < 1e-5, uniqueN(gene)])
}
br <- merge(rbindlist(br), rbindlist(h2freq), by = "cohort", all.x = TRUE)
fwrite(br, out_path(cfg, "mapt", "breadth.tsv"), sep = "\t")

print(cpgs[, .(cpg_id, cpg_pos, iv_r2 = round(iv_r2_h2_eur, 2), marg_z = round(marg_lead_z, 1), h2_z = round(h2_mqtl_z, 1),
               cond_z = round(cond_lead_z, 1), cond_n_gw, cond_mr_p = signif(cond_mr_p, 2), cond_gwas_p = signif(cond_gwas_p, 2),
               cond_H4 = round(cond_abf_H4, 2), cond_H3 = round(cond_abf_H3, 2), r2_h2 = round(r2_h2, 2), mapt_class, primary_fdr = signif(primary_fdr, 2))])
print(gw)
print(br)
print(cpgs[, .N, by = mapt_class])
