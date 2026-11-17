#!/usr/bin/env Rscript
# Phase 1, step 11: primary-hit CpGs that share one instrument.
#
# Problem. With one cis instrument per CpG, neighbouring CpGs whose instruments are the
# same SNP (or SNPs in strong LD) get Wald ratios that are the same PD signal rescaled by
# each CpG's mQTL effect. MR cannot say which CpG, if any, mediates: the exclusion
# restriction for CpG A is violated by the path through CpG B. Multivariable MR would
# need at least as many independent instruments as CpGs, which these groups don't have.
#
# What this step does, for primary hits (meta_controls x nalls2019, FDR < cfg$fdr):
#   1. Groups CpGs (same chromosome) whose lead instruments are identical or in LD
#      r^2 >= R2_GROUP (1000G EUR); single-linkage.
#   2. Per multi-CpG group:
#      - co-methylation: correlation of covariate-residualised methylation (same
#        covariates as the mQTL model) in PEG1 and PPMI controls, summarised per group;
#      - mQTL-mQTL coloc.abf for each CpG pair (controls-meta regional stats from step 5,
#        p12 = 5e-6): do the CpGs share one causal variant?
#      - MVMR feasibility: mutually independent (r^2 < 0.01) instruments in the group vs
#        group size.
#   3. Annotation (Illumina 450K manifest, hg19): gene, gene feature, CpG-island relation.
#   4. A representative CpG per group, only for reporting order: best coloc with PD
#      (SuSiE H4, else abf H4), then strongest mQTL at the instrument. It is NOT evidence
#      that this CpG is the causal one; the group is the unit of nomination.
# Summary-level output only: results/phase1/shared_iv/{cpgs,groups}.tsv
# Usage:  Rscript scripts/R/11_shared_instruments.R   (run before step 9)

source("scripts/R/phase1_lib.R")
suppressPackageStartupMessages(library(coloc))
cfg <- load_config()
R2_GROUP <- 0.8
COHORTS <- c("peg1_controls", "ppmi_controls")

hits <- fread(out_path(cfg, "mr", "meta_controls__nalls2019.tsv.gz"), colClasses = list(character = "cpg_chr"))[fdr < cfg$fdr]
inst <- fread(out_path(cfg, "instruments", "meta_controls__nalls2019.tsv.gz"),
              colClasses = list(character = c("snp_id", "chr", "cpg_chr")))[cpg_id %in% hits$cpg_id]
lead <- inst[order(-neglog10p)][!duplicated(cpg_id), .(cpg_id, cpg_chr, cpg_pos, iv_snp = snp_id, iv_pos = pos, iv_z = beta / se)]
logmsg("primary hits: %d CpGs, %d distinct lead instruments", nrow(lead), uniqueN(lead$iv_snp))

# ── 1. LD among lead instruments (1000G EUR) and grouping ────────────────────
panel <- fread(cfg$clump$panel)
keep_f <- tempfile(); writeLines(c("#IID", panel[super_pop == cfg$clump$super_pop, sample]), keep_f)
alle <- unique(inst[, .(snp_id, chr, pos, ref, alt)])
r2 <- list()
for (ch in unique(lead$cpg_chr)) {
  s <- alle[snp_id %in% inst[cpg_chr == ch, snp_id]]   # all instruments (for MVMR feasibility)
  if (nrow(s) < 2) next
  pref <- sub("\\{chrom\\}", ch, cfg$clump$ref_prefix)
  pv <- fread(cmd = sprintf("grep -v '^##' %s.pvar | cut -f1-5", pref), colClasses = list(character = c("#CHROM", "ID")))
  setnames(pv, c("chr", "pos", "kid", "kref", "kalt"))
  m <- merge(s, pv[, !"chr"], by = "pos")[(ref == kref & alt == kalt) | (ref == kalt & alt == kref)][!duplicated(snp_id)]
  if (nrow(m) < 2) next
  idf <- tempfile(); writeLines(m$kid, idf); outp <- tempfile()
  stopifnot(system(sprintf("%s --pfile %s --keep %s --extract %s --export Av --out %s --threads 4 --memory 8000 > /dev/null 2>&1",
                           cfg$clump$plink2, pref, keep_f, idf, outp)) == 0)
  tr <- fread(paste0(outp, ".traw"))
  X <- as.matrix(tr[, -(1:6)]); rownames(X) <- m$snp_id[match(tr$SNP, m$kid)]
  X <- X[apply(X, 1, sd) > 0, , drop = FALSE]
  if (nrow(X) < 2) next
  R <- cor(t(X))^2
  r2[[ch]] <- as.data.table(as.table(R))[, .(snp1 = as.character(V1), snp2 = as.character(V2), r2 = N)][snp1 < snp2]
}
r2 <- rbindlist(r2)

# single-linkage groups over CpGs: same instrument, or instruments with r^2 >= R2_GROUP
par <- setNames(lead$cpg_id, lead$cpg_id)
find <- function(x) { while (par[[x]] != x) x <- par[[x]]; x }
unite <- function(a, b) { ra <- find(a); rb <- find(b); if (ra != rb) par[[rb]] <<- ra }
for (i in seq_len(nrow(lead))) for (j in seq_len(nrow(lead))) {
  if (j <= i || lead$cpg_chr[i] != lead$cpg_chr[j]) next
  a <- lead$iv_snp[i]; b <- lead$iv_snp[j]
  rr <- if (a == b) 1 else r2[(snp1 == min(a, b) & snp2 == max(a, b)), r2]
  if (length(rr) && isTRUE(rr >= R2_GROUP)) unite(lead$cpg_id[i], lead$cpg_id[j])
}
lead[, root := vapply(cpg_id, find, "")]
lead[, group_size := .N, by = root]
setorder(lead, cpg_chr, cpg_pos)
lead[, group_id := sprintf("G%02d", match(root, unique(root)))]
logmsg("groups: %d; multi-CpG groups: %d covering %d CpGs", uniqueN(lead$group_id),
       uniqueN(lead[group_size > 1, group_id]), lead[group_size > 1, .N])

# ── 3. Annotation (450K manifest) ────────────────────────────────────────────
suppressPackageStartupMessages(library(IlluminaHumanMethylation450kanno.ilmn12.hg19))
an <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)[lead$cpg_id, c("UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_Island")])
uniq_join <- function(x) vapply(strsplit(x, ";"), function(v) paste(unique(v[nzchar(v)]), collapse = ";"), "")
lead[, `:=`(gene = uniq_join(an$UCSC_RefGene_Name), gene_feature = uniq_join(an$UCSC_RefGene_Group),
            island = an$Relation_to_Island)]

# ── 2a. Co-methylation in controls (residualised; aggregate correlations only) ──
read_header <- function(f, sep) { con <- gzcon(file(f, "rb")); on.exit(close(con)); trimws(strsplit(readLines(con, n = 1), sep)[[1]]) }
multi <- lead[group_size > 1]
cormat <- list()
for (co in COHORTS) {
  cc <- cfg$cohorts[[co]]; ddir <- pdqtl_path(cfg, cc$pipeline, "datasets", co)
  meth_f <- file.path(ddir, "methylation.csv.gz")
  geno_f <- pdqtl_path(cfg, cc$pipeline, "genotypes", paste0(co, ".geno.tsv.gz"))
  cv <- read.csv(file.path(ddir, "covariates.csv"), row.names = 1, check.names = FALSE)
  colnames(cv) <- sub("^X", "", trimws(colnames(cv)))
  idm <- read.csv(file.path(ddir, "id_map.csv"), colClasses = "character")
  al <- idm[idm[[1]] %in% read_header(meth_f, ",")[-1] & idm[[1]] %in% colnames(cv) & idm[[2]] %in% read_header(geno_f, "\t")[-1], ]
  stopifnot(nrow(al) == cc$n_samples)
  kf <- tempfile(); writeLines(multi$cpg_id, kf)
  m <- fread(cmd = sprintf("zcat %s | awk -F, 'NR==FNR{k[$1]=1; next} FNR==1 || ($1 in k)' %s -", shQuote(meth_f), kf),
             header = TRUE, colClasses = list(character = 1))
  setnames(m, read_header(meth_f, ","))
  Y <- t(as.matrix(m[, al[[1]], with = FALSE])); colnames(Y) <- m[[1]]
  C <- cbind(1, t(as.matrix(cv[, al[[1]], drop = FALSE]))); storage.mode(C) <- "double"
  Q <- qr.Q(qr(C)); Yr <- Y - Q %*% crossprod(Q, Y)
  cormat[[co]] <- cor(Yr)
}
pair_cor <- function(a, b, co) { M <- cormat[[co]]; if (a %in% rownames(M) && b %in% rownames(M)) M[a, b] else NA_real_ }

# ── 2b. mQTL-mQTL coloc on controls-meta regional statistics ─────────────────
reg <- rbindlist(lapply(COHORTS, function(co)
  fread(out_path(cfg, "regional_mqtl", paste0(co, ".tsv.gz")), colClasses = list(character = c("cpg_id", "snp_id")))[cpg_id %in% multi$cpg_id]))
ea <- rbindlist(lapply(unique(vapply(cfg$cohorts[COHORTS], `[[`, "", "pipeline")), function(pipe) {
  e <- fread(pdqtl_path(cfg, pipe, "positions/export_allele.txt"), header = FALSE, col.names = c("id", "counted"))
  e[, snp_id := sub(":[ACGT]+:[ACGT]+$", "", id)][snp_id %in% reg$snp_id][, .(pipe, snp_id, alle = sub("^[^:]+:[^:]+:", "", id))]
}))
reg[, pipe := vapply(cfg$cohorts[cohort], `[[`, "", "pipeline")]
reg <- merge(reg, ea, by = c("pipe", "snp_id"))
reg[, ok := .N == length(COHORTS) & uniqueN(alle) == 1, by = .(cpg_id, snp_id)]
rm_ <- reg[ok == TRUE][, w := 1 / se^2][, .(beta = sum(w * beta) / sum(w), se = sqrt(1 / sum(w)),
                                           eaf = sum(n * eaf) / sum(n), n = sum(n)), by = .(cpg_id, snp_id)]
mqtl_coloc <- function(a, b) {
  d <- merge(rm_[cpg_id == a], rm_[cpg_id == b], by = "snp_id", suffixes = c("1", "2"))
  if (nrow(d) < 50) return(c(NA_real_, NA_real_))
  maf <- pmin(d$eaf1, 1 - d$eaf1)
  capture.output(r <- suppressWarnings(suppressMessages(coloc.abf(
    list(beta = d$beta1, varbeta = d$se1^2, snp = d$snp_id, type = "quant", N = d$n1[1], MAF = maf),
    list(beta = d$beta2, varbeta = d$se2^2, snp = d$snp_id, type = "quant", N = d$n2[1], MAF = maf),
    p1 = 1e-4, p2 = 1e-4, p12 = 5e-6))))
  c(r$summary[["PP.H3.abf"]], r$summary[["PP.H4.abf"]])
}

prs <- list()
for (gid in unique(multi$group_id)) {
  cp <- multi[group_id == gid, cpg_id]
  for (i in seq_along(cp)) for (j in seq_along(cp)) if (j > i) {
    cc_ <- mqtl_coloc(cp[i], cp[j])
    prs[[length(prs) + 1]] <- data.table(group_id = gid, cpg1 = cp[i], cpg2 = cp[j],
                                         dist_bp = abs(lead[cpg_id == cp[i], cpg_pos] - lead[cpg_id == cp[j], cpg_pos]),
                                         same_iv = lead[cpg_id == cp[i], iv_snp] == lead[cpg_id == cp[j], iv_snp],
                                         cor_peg1 = pair_cor(cp[i], cp[j], "peg1_controls"),
                                         cor_ppmi = pair_cor(cp[i], cp[j], "ppmi_controls"),
                                         mqtl_H3 = cc_[1], mqtl_H4 = cc_[2])
  }
}
prs <- rbindlist(prs)

# ── 4. Group summaries + representative CpG ──────────────────────────────────
ab <- fread(out_path(cfg, "coloc", "meta_controls__nalls2019.tsv.gz"))[, .(cpg_id, abf_H4 = PP.H4)]
su <- fread(out_path(cfg, "coloc_susie", "meta_controls__nalls2019.tsv.gz"))[, .(cpg_id, susie_H4 = best_H4)]
lead <- merge(merge(lead, ab, by = "cpg_id", all.x = TRUE), su, by = "cpg_id", all.x = TRUE)
lead[, pd_H4 := fifelse(is.na(susie_H4), abf_H4, susie_H4)]
lead[, rank_in_group := frank(list(-fifelse(is.na(pd_H4), -1, pd_H4), -abs(iv_z)), ties.method = "first"), by = group_id]
lead[, n_instruments_cpg := inst[, .N, by = cpg_id][match(lead$cpg_id, cpg_id), N]]
# independent instruments across a group: greedy over all the group's clumped instruments,
# dropping any in LD r^2 >= 0.01 with one already kept (SNPs without LD data count once)
n_indep <- function(snps) {
  snps <- unique(snps); kept <- character(0)
  for (s in snps) {
    rr <- r2[(snp1 == s & snp2 %in% kept) | (snp2 == s & snp1 %in% kept), r2]
    if (!any(rr >= 0.01) && !s %in% kept) kept <- c(kept, s)
  }
  length(kept)
}
grp <- lead[group_size > 1, .(chr = cpg_chr[1], start = min(cpg_pos), end = max(cpg_pos), n_cpg = .N,
                              n_distinct_iv = uniqueN(inst[cpg_id %in% .SD$cpg_id, snp_id]),
                              n_independent_iv = n_indep(inst[cpg_id %in% .SD$cpg_id][order(-neglog10p), snp_id]),
                              genes = paste(unique(unlist(strsplit(gene[nzchar(gene)], ";"))), collapse = ";"),
                              representative = cpg_id[rank_in_group == 1]), by = group_id]
grp <- merge(grp, prs[, .(n_pairs = .N, min_mqtl_H4 = min(mqtl_H4, na.rm = TRUE), median_mqtl_H4 = median(mqtl_H4, na.rm = TRUE),
                          median_cor_peg1 = median(cor_peg1, na.rm = TRUE), median_cor_ppmi = median(cor_ppmi, na.rm = TRUE)),
                      by = group_id], by = "group_id")
grp[, mvmr_feasible := n_independent_iv >= n_cpg]
grp[, flag_region := flag_region(chr, start)]

fwrite(lead[, .(cpg_id, cpg_chr, cpg_pos, group_id, group_size, rank_in_group, iv_snp, iv_z, n_instruments_cpg,
                gene, gene_feature, island, pd_H4)][order(cpg_chr, cpg_pos)],
       out_path(cfg, "shared_iv", "cpgs.tsv"), sep = "\t")
fwrite(grp, out_path(cfg, "shared_iv", "groups.tsv"), sep = "\t")
fwrite(prs, out_path(cfg, "shared_iv", "pairs.tsv"), sep = "\t")
print(grp)
print(prs[, .(group_id, dist_bp, same_iv, cor_peg1 = round(cor_peg1, 2), cor_ppmi = round(cor_ppmi, 2), mqtl_H4 = round(mqtl_H4, 2))])
