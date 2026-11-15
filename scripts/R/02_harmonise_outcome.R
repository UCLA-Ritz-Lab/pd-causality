#!/usr/bin/env Rscript
# Phase 1, step 2: harmonise exposure SNPs with each PD GWAS (no rsIDs: chr:pos + alleles).
#
# Exposure effect allele = ALT, other = REF (HRC-imputed .pvar, GRCh37 forward strand).
# For each SNP, the GWAS row at the same chr:pos is matched on alleles:
#   ea/oa == ALT/REF                -> aligned          (sign +1)
#   ea/oa == REF/ALT                -> swapped          (sign -1)
#   complement of either            -> strand flip      (sign as above)
#   palindromic (A/T, C/G)          -> orientation from allele frequency; dropped if
#                                      exposure MAF >= palindrome_maf_max (ambiguous)
# After alignment, SNPs whose exposure vs GWAS ALT frequency differ by more than
# eaf_discordance_max are dropped (catches allele/strand errors; the GWAS is European).
#
# Harmonising happens before clumping, so clumping (step 3) only chooses among SNPs
# that exist in the GWAS. Availability is not an outcome-p-value criterion, so this
# does not bias the MR estimate; the number of exposure SNPs missing from the GWAS
# is reported.
#
# Output: results/phase1/harmonised/{exposure}__{outcome}.tsv.gz  (one row per SNP)
#         results/phase1/harmonised/summary.tsv
# Usage:  Rscript scripts/R/02_harmonise_outcome.R

source("scripts/R/phase1_lib.R")
cfg <- load_config()
ic  <- cfg$instruments

exp_files <- list.files(out_path(cfg, "exposure"), pattern = "\\.tsv\\.gz$", full.names = TRUE)
exp_snps <- lapply(setNames(exp_files, sub("\\.tsv\\.gz$", "", basename(exp_files))), function(f)
  unique(fread(f, select = c("snp_id", "chr", "pos", "ref", "alt", "eaf"),
               colClasses = list(character = c("snp_id", "chr")))))
all_pos <- unique(rbindlist(exp_snps)[, .(chr, pos)])
logmsg("exposure datasets: %s; unique SNP positions: %s", paste(names(exp_snps), collapse = ", "),
       format(nrow(all_pos), big.mark = ","))

harmonise <- function(x, g) {
  # x: exposure SNPs (snp_id, chr, pos, ref, alt, eaf); g: GWAS rows at those positions
  m <- merge(x, g, by = c("chr", "pos"), allow.cartesian = TRUE)
  m[, pal := is_palindromic(alt, ref)]
  m[, match := fcase(
    ea == alt & oa == ref, "aligned",
    ea == ref & oa == alt, "swapped",
    complement(ea) == alt & complement(oa) == ref, "strand_aligned",
    complement(ea) == ref & complement(oa) == alt, "strand_swapped",
    default = "mismatch")]
  m <- m[match != "mismatch"]
  # prefer the exact allele match if a position has several GWAS rows
  m[, pri := fifelse(match %in% c("aligned", "swapped"), 1L, 2L)]
  setorder(m, snp_id, pri)
  m <- m[!duplicated(snp_id)]
  m[, sgn := fifelse(match %in% c("aligned", "strand_aligned"), 1, -1)]
  # palindromic: alleles cannot tell strand; use frequency (assumes both on forward strand
  # unless frequencies say otherwise)
  m[pal == TRUE, sgn := fifelse((eaf - 0.5) * (eaf_g - 0.5) > 0, 1, -1)]
  m[, `:=`(beta_y = sgn * beta_g, eaf_y = fifelse(sgn > 0, eaf_g, 1 - eaf_g))]
  m[, maf_x := pmin(eaf, 1 - eaf)]
  m[, drop := fcase(
    pal & maf_x >= ic$palindrome_maf_max, "palindromic_ambiguous",
    abs(eaf - eaf_y) > ic$eaf_discordance_max, "eaf_discordant",
    default = NA_character_)]
  m
}

summ <- list()
for (oc in names(cfg$outcomes)) {
  o <- cfg$outcomes[[oc]]
  logmsg("=== outcome %s: reading %s", oc, o$file)
  g <- fread(o$file, select = 1:8, colClasses = list(character = "chromosome"))
  setnames(g, c("chr", "pos", "ea", "oa", "beta_g", "se_y", "eaf_g", "p_y"))
  g[, `:=`(ea = toupper(ea), oa = toupper(oa), pos = as.integer(pos))]
  n_gwas <- nrow(g)
  g <- g[all_pos, on = c("chr", "pos"), nomatch = NULL]
  logmsg("  GWAS rows: %s; at exposure positions: %s", format(n_gwas, big.mark = ","),
         format(nrow(g), big.mark = ","))

  for (ds in names(exp_snps)) {
    x <- exp_snps[[ds]]
    stopifnot(!anyDuplicated(x$snp_id))
    m <- harmonise(x, g)
    at_pos <- x[g, on = c("chr", "pos"), nomatch = NULL, unique(snp_id)]
    keep <- m[is.na(drop)]
    fwrite(keep[, .(snp_id, chr, pos, ref, alt, eaf_x = eaf, beta_y, se_y, eaf_y, p_y,
                    harmonise = match, palindromic = pal)],
           out_path(cfg, "harmonised", sprintf("%s__%s.tsv.gz", ds, oc)), sep = "\t")
    s <- data.table(exposure = ds, outcome = oc, snps_exposure = nrow(x),
                    at_gwas_position = length(at_pos), allele_matched = nrow(m),
                    allele_mismatch = length(at_pos) - nrow(m),
                    strand_flipped = m[match %like% "strand", .N],
                    palindromic_kept = keep[pal == TRUE, .N],
                    dropped_palindromic = m[drop %in% "palindromic_ambiguous", .N],
                    dropped_eaf_discordant = m[drop %in% "eaf_discordant", .N],
                    harmonised = nrow(keep))
    s[, pct_missing_from_gwas := round(100 * (1 - harmonised / snps_exposure), 1)]
    summ[[paste(ds, oc)]] <- s
    logmsg("  %s: %s of %s SNPs harmonised (%.1f%% missing/dropped)", ds,
           format(nrow(keep), big.mark = ","), format(nrow(x), big.mark = ","), s$pct_missing_from_gwas)
  }
  rm(g); gc(verbose = FALSE)
}
s <- rbindlist(summ)
fwrite(s, out_path(cfg, "harmonised", "summary.tsv"), sep = "\t")
print(s)
