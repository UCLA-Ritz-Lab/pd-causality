#!/usr/bin/env Rscript
# Phase 1, step 8: in-house PD GWAS in PPMI and PEG for the crossed two-sample design
# (PEG mQTL x PPMI GWAS, PPMI mQTL x PEG GWAS; never same-cohort pairs).
#
# PLINK2 Firth-penalised logistic regression for every variant on imputed dosages (HRC r1.1, hg19),
# variants with MAF >= 0.01 (in the analysed sample) and imputation R2 >= 0.3.
# Firth for all variants, not just as a fallback: plink2 falls back only when the ordinary
# fit fails to converge, so quasi-separated rare variants (e.g. LRRK2 G2019S, carried almost
# only by cases) came out with |log-OR| in the hundreds. Firth gives finite, shrunk estimates.
# European-ancestry participants only (ANCESTRY_REF assigned_superpop == EUR), to match
# Nalls and the mQTL instruments. Covariates: sex, age, ancestry PCs 1-5.
#   ppmi       : PD (all 511 genotyped PD) vs healthy controls; prodromal + SWEDD excluded.
#                age = ENROLL_AGE; sex from Demographics (PPMI coding 1 = male).
#   ppmi_nogen : as ppmi but without PD cases enrolled through LRRK2/GBA/SNCA
#                (ascertained on genotype).
#   peg        : PEG1 + PEG2 participants with covariates (the methylation study; the
#                full PEG genotype set has no phenotype file locally). PD status
#                PDstudyDiseaseNumeric; age, Female from the PEG covariate files.
# Individual-level phenotype/covariate files are written only under results/ (git-ignored).
#
# Output: results/phase1/inhouse_gwas/{name}.tsv.gz in Nalls column format
#   (chromosome, base_pair_location, effect_allele, other_allele, beta = log-OR,
#    standard_error, effect_allele_frequency, p_value, n_cases, n_controls)
# Usage:  Rscript scripts/R/08_inhouse_gwas.R [name ...]

source("scripts/R/phase1_lib.R")
cfg  <- load_config()
args <- commandArgs(trailingOnly = TRUE)
odir <- out_path(cfg, "inhouse_gwas")
dir.create(odir, showWarnings = FALSE)
anc <- fread(pdqtl_path(cfg, "ANCESTRY_REF/ancestry_pcs_all.csv"), select = c("IID", "cohort", "assigned_superpop"))

build_ppmi <- function(exclude_genetic) {
  ps <- fread("data/raw/PPMI_covariates/Participant_Status_14Dec2025.csv", colClasses = list(character = "PATNO"))
  dm <- fread("data/raw/PPMI_covariates/Demographics_14Dec2025.csv", select = c("PATNO", "SEX"), colClasses = list(character = "PATNO"))
  pcs <- fread(pdqtl_path(cfg, "ANCESTRY_REF/covariates/ancestry_pcs_PPMI.csv"))
  d <- ps[COHORT_DEFINITION %in% c("Parkinson's Disease", "Healthy Control")]
  if (exclude_genetic) d <- d[!(COHORT_DEFINITION == "Parkinson's Disease" & (ENRLLRRK2 == 1 | ENRLGBA == 1 | ENRLSNCA == 1))]
  d <- merge(d, dm, by = "PATNO")
  d[, `:=`(IID = paste0("PPMISI", PATNO), PD = fifelse(COHORT_DEFINITION == "Parkinson's Disease", 2L, 1L),
           MALE = SEX, AGE = ENROLL_AGE)]
  d <- merge(d[, .(IID, PD, MALE, AGE)], pcs[, .(IID, PC1, PC2, PC3, PC4, PC5)], by = "IID")
  d[IID %in% anc[cohort == "PPMI" & assigned_superpop == "EUR", IID]]
}

build_peg <- function() {
  lk <- fread("data/raw/PEG_WGS/gwas_linkage.csv", colClasses = "character")[DUP == "0"]
  c1 <- fread("data/raw/PEG_covariates/covariates_peg1.csv", colClasses = list(character = "ExternalDNACode"))
  c2 <- fread("data/raw/PEG_covariates/covariates_peg2.csv", colClasses = list(character = "PEGID"))
  cv <- unique(rbind(c1[, .(Pegid = ExternalDNACode, PDn = PDstudyDiseaseNumeric, Female, Age)],
                     c2[, .(Pegid = PEGID, PDn = PDstudyDiseaseNumeric, Female, Age)]), by = "Pegid")
  pcs <- fread(pdqtl_path(cfg, "ANCESTRY_REF/covariates/ancestry_pcs_PEG.csv"))
  d <- merge(cv, unique(lk[, .(Pegid, IID = GWAS_ID)], by = "Pegid"), by = "Pegid")
  d[, `:=`(PD = PDn + 1L, MALE = 1L - Female, AGE = Age)]
  d <- merge(d[, .(IID, PD, MALE, AGE)], pcs[, .(IID, PC1, PC2, PC3, PC4, PC5)], by = "IID")
  d <- d[!duplicated(IID)]
  d[IID %in% anc[cohort == "PEG" & assigned_superpop == "EUR", IID]]
}

runs <- list(
  ppmi       = list(pfile = pdqtl_path(cfg, "PPMI_WGS/post_imputed/ppmi_imputed"), build = function() build_ppmi(FALSE)),
  ppmi_nogen = list(pfile = pdqtl_path(cfg, "PPMI_WGS/post_imputed/ppmi_imputed"), build = function() build_ppmi(TRUE)),
  peg        = list(pfile = pdqtl_path(cfg, "PEG_WGS/post_imputed/peg_imputed"),  build = build_peg))
if (length(args)) runs <- runs[args]

for (nm in names(runs)) {
  r <- runs[[nm]]
  d <- r$build()
  psam <- fread(paste0(r$pfile, ".psam"), colClasses = "character")
  d <- d[IID %in% psam$IID & complete.cases(d)]
  logmsg("%s: %d cases, %d controls (EUR, complete covariates)", nm, d[PD == 2, .N], d[PD == 1, .N])
  pf <- file.path(odir, paste0(nm, ".pheno"))
  fwrite(d[, .(`#FID` = IID, IID, PD)], pf, sep = "\t")   # psam FID == IID
  fwrite(d[, .(`#FID` = IID, IID, MALE, AGE, PC1, PC2, PC3, PC4, PC5)], file.path(odir, paste0(nm, ".covar")), sep = "\t")
  outp <- file.path(odir, nm)
  cmd <- sprintf(paste("%s --pfile %s --keep %s --pheno %s --pheno-name PD --covar %s --covar-variance-standardize",
                       "--extract-if-info 'R2 >= 0.3' --maf 0.01 --glm firth hide-covar cols=+a1freq",
                       "--threads 4 --memory 8000 --out %s > %s.plink.log 2>&1"),
                 cfg$clump$plink2, r$pfile, pf, pf, file.path(odir, paste0(nm, ".covar")), outp, outp)
  gpat <- paste0("^", nm, "\\.PD\\.glm\\.(logistic|logistic\\.hybrid|firth)$")
  if (!length(list.files(odir, pattern = gpat))) {        # reuse a finished plink2 run
    t0 <- Sys.time()
    stopifnot(system(cmd) == 0)
    logmsg("  plink2 done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
  }
  gf <- list.files(odir, pattern = gpat, full.names = TRUE)
  stopifnot(length(gf) == 1)
  g <- fread(gf, colClasses = list(character = "#CHROM"))
  setnames(g, "#CHROM", "CHROM")
  for (cn in c("OR", "LOG(OR)_SE", "P", "A1_FREQ")) set(g, j = cn, value = suppressWarnings(as.numeric(g[[cn]])))
  n_fail <- g[TEST == "ADD" & !is.finite(`LOG(OR)_SE`), .N]
  g <- g[TEST == "ADD" & is.finite(`LOG(OR)_SE`) & is.finite(OR) & OR > 0]
  logmsg("  dropped %s variants with failed fits (NA)", format(n_fail, big.mark = ","))
  # A1 is the tested allele; other allele = whichever of REF/ALT is not A1
  g[, other := fifelse(A1 == ALT, REF, ALT)]
  out <- g[, .(chromosome = CHROM, base_pair_location = POS, effect_allele = A1, other_allele = other,
               beta = log(OR), standard_error = `LOG(OR)_SE`, effect_allele_frequency = A1_FREQ,
               p_value = P, n_cases = d[PD == 2, .N], n_controls = d[PD == 1, .N])]
  fwrite(out, file.path(odir, paste0(nm, ".tsv.gz")), sep = "\t")
  unlink(gf)
  lam <- median(qchisq(out$p_value, 1, lower.tail = FALSE), na.rm = TRUE) / qchisq(0.5, 1)
  logmsg("  %s variants written; lambda_GC = %.3f", format(nrow(out), big.mark = ","), lam)
}
