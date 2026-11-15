#!/usr/bin/env Rscript
# Phase 1, step 3: LD-clump candidate instruments separately for each CpG.
#
# Clumping must be per CpG (the same SNP can instrument several CpGs), so instead of
# calling plink --clump thousands of times we export 1000G EUR dosages for every
# harmonised candidate SNP once per chromosome and clump greedily in R:
#   repeat: index = most significant remaining SNP (exposure p) present in the LD
#           reference; drop remaining SNPs with r^2 >= clump$r2 with the index.
# r^2 is computed among all candidates of the CpG (cis window 1 Mb), with no distance cap.
#
# LD reference: 1000G phase 3 GRCh37 EUR (n = 503), from the ANCESTRY_REF panel
# (biallelic SNPs with global MAF >= 5%). SNPs absent from it, or monomorphic in EUR,
# cannot be clumped and are not eligible as index SNPs. A CpG with no eligible SNP keeps
# its single most significant SNP (flag ld_ref = "missing").
#
# Output: results/phase1/instruments/{exposure}__{outcome}.tsv.gz  (instrument rows,
#         exposure + outcome effects aligned to ALT), and summary.tsv
# Usage:  Rscript scripts/R/03_clump.R

source("scripts/R/phase1_lib.R")
cfg <- load_config()
cc  <- cfg$clump
tmp <- out_path(cfg, "tmp_ld")
dir.create(tmp, showWarnings = FALSE)

# ── Combos and candidates ────────────────────────────────────────────────────
hfiles <- list.files(out_path(cfg, "harmonised"), pattern = "__.*\\.tsv\\.gz$", full.names = TRUE)
combos <- sub("\\.tsv\\.gz$", "", basename(hfiles))
cand <- list()
for (cb in combos) {
  ds <- sub("__.*", "", cb)
  h  <- fread(file.path(out_path(cfg, "harmonised"), paste0(cb, ".tsv.gz")),
              colClasses = list(character = c("snp_id", "chr")))
  e  <- fread(out_path(cfg, "exposure", paste0(ds, ".tsv.gz")),
              colClasses = list(character = c("snp_id", "chr", "cpg_chr")))
  cand[[cb]] <- merge(e, h[, !c("chr", "pos", "ref", "alt")], by = "snp_id")
  logmsg("%s: %s candidate pairs, %s CpGs", cb, format(nrow(cand[[cb]]), big.mark = ","),
         format(uniqueN(cand[[cb]]$cpg_id), big.mark = ","))
}
all_snps <- unique(rbindlist(lapply(cand, function(d) d[, .(snp_id, chr, pos, ref, alt)])))

# ── EUR keep list ────────────────────────────────────────────────────────────
panel <- fread(cc$panel)
keep_f <- file.path(tmp, "eur.keep")
writeLines(c("#IID", panel[super_pop == cc$super_pop, sample]), keep_f)

# ── Clump one CpG ────────────────────────────────────────────────────────────
clump_cpg <- function(ids, G) {
  # ids: SNP ids ordered by significance; G: standardised dosage rows (or NULL rows absent)
  ok <- ids[ids %in% rownames(G)]
  if (!length(ok)) return(list(index = ids[1], ld_ref = "missing"))
  idx <- character(0); rem <- ok
  while (length(rem)) {
    i <- rem[1]; idx <- c(idx, i)
    rem <- rem[-1]
    if (!length(rem)) break
    r <- as.vector(G[rem, , drop = FALSE] %*% G[i, ]) / (ncol(G) - 1)
    rem <- rem[r^2 < cc$r2]
  }
  list(index = idx, ld_ref = if (ok[1] == ids[1]) "lead_in_ref" else "lead_not_in_ref")
}

res <- setNames(vector("list", length(combos)), combos)
for (ch in c(as.character(1:22))) {
  s <- all_snps[chr == ch]
  if (!nrow(s)) next
  pref <- sub("\\{chrom\\}", ch, cc$ref_prefix)
  pvar <- fread(cmd = sprintf("grep -v '^##' %s.pvar | cut -f1-5", pref),
                colClasses = list(character = c("#CHROM", "ID")))
  setnames(pvar, c("chr", "pos", "ref_id", "kref", "kalt"))
  # match on position + same allele pair (either orientation; r^2 ignores sign)
  m <- merge(s, pvar, by = "pos")[(ref == kref & alt == kalt) | (ref == kalt & alt == kref)]
  m <- m[!duplicated(snp_id)]
  ids_f <- file.path(tmp, sprintf("chr%s.ids", ch))
  writeLines(m$ref_id, ids_f)
  outp <- file.path(tmp, sprintf("chr%s", ch))
  cmd <- sprintf("%s --pfile %s --keep %s --extract %s --export Av --out %s --threads 4 --memory 8000 > %s.plink.log 2>&1",
                 cc$plink2, pref, keep_f, ids_f, outp, outp)
  stopifnot(system(cmd) == 0)
  tr <- fread(paste0(outp, ".traw"))
  G <- as.matrix(tr[, -(1:6)])
  rownames(G) <- m$snp_id[match(tr$SNP, m$ref_id)]
  stopifnot(!anyNA(rownames(G)))
  sds <- apply(G, 1, sd, na.rm = TRUE)
  G <- G[!is.na(sds) & sds > 0, , drop = FALSE]
  G <- (G - rowMeans(G, na.rm = TRUE)) / apply(G, 1, sd, na.rm = TRUE)
  G[is.na(G)] <- 0
  unlink(paste0(outp, ".traw"))
  logmsg("chr%s: %s candidate SNPs, %s in EUR ref (polymorphic)", ch,
         format(nrow(s), big.mark = ","), format(nrow(G), big.mark = ","))

  for (cb in combos) {
    d <- cand[[cb]][cpg_chr == ch]
    if (!nrow(d)) next
    setorder(d, cpg_id, -neglog10p, -F_stat)
    out <- d[, {
      cl <- clump_cpg(snp_id, G)
      .(snp_id = cl$index, ld_ref = cl$ld_ref, n_candidates = .N,
        n_candidates_in_ref = sum(snp_id %in% rownames(G)))
    }, by = cpg_id]
    res[[cb]][[ch]] <- merge(out, d, by = c("cpg_id", "snp_id"))
  }
  rm(G, tr); gc(verbose = FALSE)
}

summ <- list()
for (cb in combos) {
  r <- rbindlist(res[[cb]])
  r[, n_instruments := .N, by = cpg_id]
  setorder(r, cpg_id, -neglog10p)
  fwrite(r, out_path(cfg, "instruments", paste0(cb, ".tsv.gz")), sep = "\t")
  per_cpg <- unique(r[, .(cpg_id, n_instruments, ld_ref)])
  summ[[cb]] <- data.table(
    combo = cb, cpgs = nrow(per_cpg), instruments = nrow(r),
    cpgs_1_iv = per_cpg[n_instruments == 1, .N], cpgs_2_iv = per_cpg[n_instruments == 2, .N],
    cpgs_3plus_iv = per_cpg[n_instruments >= 3, .N],
    lead_in_ref = per_cpg[ld_ref == "lead_in_ref", .N],
    lead_not_in_ref = per_cpg[ld_ref == "lead_not_in_ref", .N],
    no_snp_in_ref = per_cpg[ld_ref == "missing", .N])
}
s <- rbindlist(summ)
fwrite(s, out_path(cfg, "instruments", "summary.tsv"), sep = "\t")
print(s)
