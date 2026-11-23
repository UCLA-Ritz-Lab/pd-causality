#!/usr/bin/env Rscript
# Phase 1, step 12: SNPs inside the array probe body, for the primary MR hits.
#
# Why. Step 1 drops instruments within 50 bp of the CpG, so no instrument lies in a probe.
# But an instrument can be in LD with a SNP that does. A probe-body SNP mismatches the
# 50-mer in carriers and can change the measured beta-value without any change in
# methylation (hybridisation artefact). The mQTL is then partly an array artefact, and
# because the two SNPs are in LD, nothing in the array data can separate the two.
#
# Probe body. The 50-mer (ProbeSeqA from minfi's 450K manifest, EPIC manifest for
# EPIC-only probes) is matched to the bisulfite-converted hg19 reference (+-80 bp around
# the CpG, both strands, methylated and unmethylated conversion, <= 3 mismatches). Every
# probe must match exactly once. For Type I probes the single-base-extension (SBE) base
# next to the probe's 3' end is reported too.
#
# SNPs. Candidate SNPs inside the probe body: any SNP in the step 5 regional mQTL files
# (all cohorts) or in the 1000G EUR reference. For each, EUR r^2 with every instrument of
# the CpG (1000G phase 3 EUR; SNPs absent from the reference get NA) and its strongest
# regional mQTL p in the primary-exposure cohorts.
#
# Flag (reported, does not change the tier): probe_snp_ld = TRUE if any probe-body SNP
# has r^2 >= 0.8 with an instrument of the CpG.
#
# Output: results/phase1/probe_snps/cpgs.tsv (one row per primary hit),
#         results/phase1/probe_snps/snps.tsv (one row per CpG x probe-body SNP).
# Summary-level only. Usage: Rscript scripts/R/12_probe_snps.R

source("scripts/R/phase1_lib.R")
suppressPackageStartupMessages({
  library(Biostrings); library(Rsamtools); library(GenomicRanges); library(minfi)
  library(IlluminaHumanMethylation450kmanifest); library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  library(IlluminaHumanMethylationEPICmanifest); library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
})
cfg <- load_config()
R2_FLAG <- 0.8
WIN <- 80
fa <- FaFile(cfg$probe_check$fasta)
fa_seqs <- seqlevels(scanFaIndex(fa))

hits <- fread(out_path(cfg, "mr", "meta_controls__nalls2019.tsv.gz"), colClasses = list(character = "cpg_chr"))[fdr < cfg$fdr]
inst <- fread(out_path(cfg, "instruments", "meta_controls__nalls2019.tsv.gz"),
              colClasses = list(character = c("snp_id", "chr", "cpg_chr")))[cpg_id %in% hits$cpg_id]
logmsg("primary hits: %d CpGs, %d instruments", nrow(hits), nrow(inst))

# ── probe sequences and positions (450K first, EPIC for EPIC-only probes) ────
probe_tab <- function(man, ann) {
  pr <- rbind(as.data.table(as.data.frame(getProbeInfo(man, type = "I")))[, .(cpg_id = Name, type = "I", seq = ProbeSeqA)],
              as.data.table(as.data.frame(getProbeInfo(man, type = "II")))[, .(cpg_id = Name, type = "II", seq = ProbeSeqA)])
  loc <- as.data.table(as.data.frame(getLocations(ann)), keep.rownames = "cpg_id")[, .(cpg_id, chr = as.character(seqnames), pos = start)]
  merge(pr, loc, by = "cpg_id")
}
pt <- probe_tab(IlluminaHumanMethylation450kmanifest, IlluminaHumanMethylation450kanno.ilmn12.hg19)[cpg_id %in% hits$cpg_id][, array := "450K"]
miss <- setdiff(hits$cpg_id, pt$cpg_id)
if (length(miss))
  pt <- rbind(pt, probe_tab(IlluminaHumanMethylationEPICmanifest, IlluminaHumanMethylationEPICanno.ilm10b4.hg19)[cpg_id %in% miss][, array := "EPIC"])
stopifnot(all(hits$cpg_id %in% pt$cpg_id))

bisulfite <- function(s, meth) {
  s <- as.character(s)
  if (meth) { s <- gsub("CG", "XG", s); s <- gsub("C", "T", s); gsub("X", "C", s) } else gsub("C", "T", s)
}
map_probe <- function(chr, pos, seq) {
  ch <- if (chr %in% fa_seqs) chr else sub("^chr", "", chr)
  g0 <- pos - WIN
  gs <- getSeq(fa, GRanges(ch, IRanges(g0, pos + WIN)))[[1]]
  L <- length(gs); ps <- DNAString(seq); hit <- NULL
  for (tm in c("+", "-")) for (meth in c(TRUE, FALSE)) {
    conv <- DNAString(bisulfite(if (tm == "+") gs else reverseComplement(gs), meth))
    for (pp in list(ps, reverseComplement(ps))) {
      m <- matchPattern(pp, conv, max.mismatch = 3, fixed = FALSE)
      for (k in seq_along(m)) {
        st <- start(m)[k]; en <- end(m)[k]
        if (tm == "-") { a <- L - en + 1; en <- L - st + 1; st <- a }
        hit <- rbind(hit, data.table(st = g0 + st - 1, en = g0 + en - 1))
      }
    }
  }
  hit <- unique(hit)
  if (is.null(hit) || nrow(hit) != 1) stop("probe at ", chr, ":", pos, " matched ", if (is.null(hit)) 0 else nrow(hit), " times")
  hit
}
pt[, c("probe_start", "probe_end") := {
  h <- map_probe(chr, pos, seq); list(h$st, h$en)
}, by = cpg_id]
# Type I: SBE base is the one just beyond the probe's 3' end, which sits at the CpG
pt[, sbe_pos := fifelse(type == "I", fifelse(probe_start >= pos - 1, probe_start - 1L, probe_end + 1L), NA_integer_)]
pt[, cpg_chr := sub("^chr", "", chr)]
logmsg("probes mapped: %d (450K %d, EPIC %d)", nrow(pt), pt[array == "450K", .N], pt[array == "EPIC", .N])

# ── candidate SNPs in the probe body ─────────────────────────────────────────
reg <- rbindlist(lapply(names(cfg$cohorts), function(co) {
  f <- out_path(cfg, "regional_mqtl", paste0(co, ".tsv.gz"))
  if (!file.exists(f)) return(NULL)
  fread(f, colClasses = list(character = c("cpg_id", "snp_id")))[cpg_id %in% hits$cpg_id & abs(dist_bp) <= 60,
                                                                 .(cpg_id, snp_id, dist_bp, p, cohort)]
}))
prim_cohorts <- unlist(cfg$exposures$meta_controls$cohorts)
reg_best <- reg[cohort %in% prim_cohorts, .(mqtl_p_min = min(p)), by = .(cpg_id, snp_id)]

panel <- fread(cfg$clump$panel)
keep_f <- tempfile(); writeLines(c("#IID", panel[super_pop == cfg$clump$super_pop, sample]), keep_f)

snps <- list()
for (ch in unique(pt$cpg_chr)) {
  pref <- sub("\\{chrom\\}", ch, cfg$clump$ref_prefix)
  pv <- fread(cmd = sprintf("grep -v '^##' %s.pvar | cut -f1-5", pref), colClasses = list(character = c("#CHROM", "ID")))
  setnames(pv, c("chr", "pos", "kid", "kref", "kalt"))
  pc <- pt[cpg_chr == ch]
  cand <- rbindlist(lapply(seq_len(nrow(pc)), function(i) {
    r <- pc[i]
    a <- reg[cpg_id == r$cpg_id, .(snp_id, pos = as.integer(sub(".*:", "", snp_id)))]
    b <- pv[pos >= r$probe_start & pos <= r$probe_end, .(snp_id = paste0(ch, ":", pos), pos)]
    x <- unique(rbind(a, b))[pos >= r$probe_start & pos <= r$probe_end]
    if (!nrow(x)) return(NULL)
    x[, `:=`(cpg_id = r$cpg_id, cpg_pos = r$pos, sbe_pos = r$sbe_pos, probe_type = r$type)]
  }))
  if (!nrow(cand)) next
  cand[, position := fcase(pos %in% c(cpg_pos, cpg_pos + 1L), "CpG",
                           !is.na(sbe_pos) & pos == sbe_pos, "SBE", default = "body")]
  # EUR r^2 between each probe-body SNP and the CpG's instruments
  iv <- inst[cpg_id %in% cand$cpg_id, .(cpg_id, iv_snp = snp_id, iv_pos = pos)]
  ids <- unique(pv[pos %in% unique(c(cand$pos, iv$iv_pos))][!duplicated(pos)])
  dos <- NULL
  if (nrow(ids) >= 2) {
    idf <- tempfile(); writeLines(ids$kid, idf); outp <- tempfile()
    stopifnot(system(sprintf("%s --pfile %s --keep %s --extract %s --export Av --out %s --threads 4 --memory 8000 > /dev/null 2>&1",
                             cfg$clump$plink2, pref, keep_f, idf, outp)) == 0)
    tr <- fread(paste0(outp, ".traw"))
    dos <- as.matrix(tr[, -(1:6)]); rownames(dos) <- tr$POS
  }
  r2_of <- function(a, b) {
    a <- as.character(a); b <- as.character(b)
    if (is.null(dos) || !a %in% rownames(dos) || !b %in% rownames(dos)) return(NA_real_)
    if (a == b) return(1)
    x <- dos[a, ]; y <- dos[b, ]
    if (sd(x) == 0 || sd(y) == 0) return(NA_real_)
    cor(x, y)^2
  }
  cand <- merge(cand, iv, by = "cpg_id", allow.cartesian = TRUE)
  cand[, r2_iv := mapply(r2_of, pos, iv_pos)]
  cand <- cand[order(-r2_iv)][, .(dist_bp = pos[1] - cpg_pos[1], position = position[1], probe_type = probe_type[1],
                                  max_r2_iv = suppressWarnings(max(r2_iv, na.rm = TRUE)),
                                  iv_best = iv_snp[1], in_1000g = pos[1] %in% ids$pos), by = .(cpg_id, snp_id)]
  cand[!is.finite(max_r2_iv), max_r2_iv := NA_real_]
  snps[[ch]] <- cand
}
snps <- rbindlist(snps)
snps <- merge(snps, reg_best, by = c("cpg_id", "snp_id"), all.x = TRUE)
# regional lead (primary-exposure cohorts, strongest p) for each CpG
lead <- rbindlist(lapply(prim_cohorts, function(co)
  fread(out_path(cfg, "regional_mqtl", paste0(co, ".tsv.gz")), colClasses = list(character = c("cpg_id", "snp_id")))[cpg_id %in% hits$cpg_id]))
lead <- lead[order(p)][!duplicated(cpg_id), .(cpg_id, regional_lead = snp_id)]
snps <- merge(snps, lead, by = "cpg_id", all.x = TRUE)
snps[, is_regional_lead := snp_id == regional_lead][, regional_lead := NULL]
setorder(snps, cpg_id, dist_bp)
fwrite(snps, out_path(cfg, "probe_snps", "snps.tsv"), sep = "\t")

# ── per-CpG summary ──────────────────────────────────────────────────────────
cs <- snps[, .(n_probe_snps = .N, probe_snps = paste(snp_id, collapse = ";"),
               probe_snp_positions = paste(unique(position), collapse = ";"),
               probe_snp_max_r2_iv = suppressWarnings(max(max_r2_iv, na.rm = TRUE)),
               probe_snp_min_mqtl_p = suppressWarnings(min(mqtl_p_min, na.rm = TRUE)),
               probe_snp_is_regional_lead = any(is_regional_lead)), by = cpg_id]
cs[!is.finite(probe_snp_max_r2_iv), probe_snp_max_r2_iv := NA_real_]
cs[!is.finite(probe_snp_min_mqtl_p), probe_snp_min_mqtl_p := NA_real_]
out <- merge(pt[, .(cpg_id, array, probe_type = type, probe_start, probe_end, sbe_pos)], cs, by = "cpg_id", all.x = TRUE)
out[is.na(n_probe_snps), `:=`(n_probe_snps = 0L, probe_snp_is_regional_lead = FALSE)]
out[, probe_snp_ld := !is.na(probe_snp_max_r2_iv) & probe_snp_max_r2_iv >= R2_FLAG]
fwrite(out, out_path(cfg, "probe_snps", "cpgs.tsv"), sep = "\t")
logmsg("CpGs with >= 1 probe-body SNP: %d; probe_snp_ld (r2 >= %.1f with an instrument): %d; probe SNP is regional lead: %d",
       out[n_probe_snps > 0, .N], R2_FLAG, out[probe_snp_ld == TRUE, .N], out[probe_snp_is_regional_lead == TRUE, .N])
