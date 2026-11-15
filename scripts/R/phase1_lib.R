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
