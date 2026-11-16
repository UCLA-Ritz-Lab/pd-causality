# Shared helpers for Phase 1 (two-sample MR: cis-mQTL exposure x PD GWAS outcome).
# Source from the repo root: source("scripts/R/phase1_lib.R")

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})

load_config <- function(path = "config/config.yaml") {
  if (!file.exists(path)) stop("Missing ", path, " - copy config/config.example.yaml and edit it.")
  cfg <- yaml::read_yaml(path)
  dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
  cfg
}

pdqtl_path <- function(cfg, ...) file.path(cfg$pdqtl, ...)
out_path   <- function(cfg, ...) { p <- file.path(cfg$out_dir, ...); dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE); p }

logmsg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%F %T"), sprintf(...)))

# -log10 two-sided p from z without underflow (p < 1e-308 is common for cis-mQTLs).
neglog10p_z <- function(z) -(pnorm(-abs(z), log.p = TRUE) + log(2)) / log(10)

# z threshold equivalent to a two-sided p threshold.
z_for_p <- function(p) qnorm(p / 2, lower.tail = FALSE)

complement <- function(x) chartr("ACGT", "TGCA", x)
is_palindromic <- function(a1, a2) a1 == complement(a2)

# Nalls-format GWAS (first 8 columns), alleles uppercased (GCST009324 is lowercase).
read_gwas <- function(f) {
  g <- fread(f, select = 1:8, colClasses = list(character = "chromosome"))
  setnames(g, c("chr", "pos", "ea", "oa", "beta_g", "se_y", "eaf_g", "p_y"))
  g[, `:=`(ea = toupper(ea), oa = toupper(oa), pos = as.integer(pos))]
  g
}

# Align GWAS rows to the exposure's ALT allele.
# x: exposure SNPs (snp_id, chr, pos, ref, alt, eaf = ALT freq); g: read_gwas() rows.
# ic: config$instruments (palindrome_maf_max, eaf_discordance_max).
# Returns one row per matched snp_id with beta_y / eaf_y on the ALT allele and `drop`
# set (non-NA) for ambiguous palindromes and allele-frequency-discordant SNPs.
harmonise_gwas <- function(x, g, ic) {
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
  # palindromic: alleles cannot tell strand; orient by allele frequency
  m[pal == TRUE, sgn := fifelse((eaf - 0.5) * (eaf_g - 0.5) > 0, 1, -1)]
  m[, `:=`(beta_y = sgn * beta_g, eaf_y = fifelse(sgn > 0, eaf_g, 1 - eaf_g))]
  m[, maf_x := pmin(eaf, 1 - eaf)]
  m[, drop := fcase(
    pal & maf_x >= ic$palindrome_maf_max, "palindromic_ambiguous",
    abs(eaf - eaf_y) > ic$eaf_discordance_max, "eaf_discordant",
    default = NA_character_)]
  m
}

# Long-range LD / pleiotropy-prone regions to flag (GRCh37). Generous windows.
# MAPT H1/H2 inversion, HLA, GBA1, LRRK2 (CLAUDE.md: treat with special care).
FLAG_REGIONS <- data.table(
  region = c("MAPT_17q21.31", "HLA", "GBA1", "LRRK2"),
  chr    = c("17", "6", "1", "12"),
  start  = c(43400000, 28477797, 154200000, 39600000),
  end    = c(44900000, 33448354, 156200000, 41800000)
)
flag_region <- function(chr, pos) {
  out <- rep(NA_character_, length(chr))
  for (i in seq_len(nrow(FLAG_REGIONS))) {
    r <- FLAG_REGIONS[i]
    out[chr == r$chr & pos >= r$start & pos <= r$end] <- r$region
  }
  out
}
