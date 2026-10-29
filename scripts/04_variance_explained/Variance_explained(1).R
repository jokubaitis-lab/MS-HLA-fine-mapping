if (!nzchar(Sys.getenv("HLA_PROJECT_ROOT", unset = ""))) stop("Set HLA_PROJECT_ROOT to the input-data root before running this historical script.")
# =============================================================================
# Variance explained by amino-acid positions (Nagelkerke pseudo-R²)
# - Per-position R² (covariates + that position's residue dummies)
# - Cumulative R² in a chosen order (delta R² at each step)
# - Partial/unique R² by drop-one from the full model
#
# Assumes you already loaded:
#   - dose: data.frame with IID and HLA allele dosage columns (0..2)
#   - allelenames: character vector of allele column names (matching dose)
#   - pheno: data.frame with IID and the binary phenotype column
#   - cov:   data.frame with IID and PCs (PC1..PC5)
#   - pheno_name: string giving the phenotype column name in `pheno`
#   - info:  the full AA dictionary (NOT filtered to one gene):
#            readRDS("HLA_DICTIONARY_AA.hg19.imgt3320.AA_tf.in_ref.rds")
# It reuses the same covariate term string and PCs as your script.
#
# If you’re running from a clean session, set RUN_MINIMAL_LOADER <- TRUE
# and edit the file paths under the "Minimal loader" block.
# =============================================================================

# -------------------------- Packages (install if missing) ---------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(tidyr)
})

options(stringsAsFactors = FALSE)
`%&%` <- function(a,b) paste0(a,b)

# -------------------------- Toggle: minimal loader ----------------------------
RUN_MINIMAL_LOADER <- TRUE  # <<< set TRUE if starting from a clean session

if (RUN_MINIMAL_LOADER) {
  # ----------- Minimal loader (EDIT these paths/placeholders) -----------------
  data_assoc  <- file.path(Sys.getenv("HLA_PROJECT_ROOT", unset = ""), "project/HLA_Imputation/HLA_analyses_tutorial/data_assoc/")
  raw_prefix  <- "AU_HLA_imputed_merged.keep.qc_filtered.QCed_HLA_tf.imputed"       # same prefix used for .raw
  imputed     <- "AU_HLA_imputed_merged.keep.qc_filtered.QCed_HLA_tf"      # same prefix used for *.allele_name.txt
  pheno_file  <- "Case_control.txt"
  cov_name    <- "Case_Control_covariates.txt"
  pheno_name  <- "MS_pheno"                       # set to your phenotype column name

  info <- readRDS(data_assoc %&% "HLA_DICTIONARY_AA.hg19.imgt3320.AA_tf.in_ref.rds")

  dose <- read.table(raw_prefix %&% ".raw", header = TRUE, check.names = FALSE)
  dose <- dose[, c(-1, -3:-6)]  # drop FID,PAT,MAT,SEX,PHENOTYPE
  allelenames <- as.character(read.table(imputed %&% ".allele_name.txt")$V1)
  colnames(dose) <- c("IID", allelenames)

  pheno <- read.table(pheno_file, header = TRUE)
  cov   <- read.table(cov_name, header = TRUE)
  cov   <- cov[, (colnames(cov) %in% c("IID","PC1","PC2","PC3","PC4","PC5","D_MSGP1_AU","D_MSGP2_AU","D_MSPGX_AU"))]
}

# -------------------------- Config: positions to evaluate ---------------------
pos_spec <- list(
  DRB1 = c("71","13","86"),
  A    = c("95"),
  DPB1 = c("65")
)

# Use the same covariate term string as your stepwise code.
if (!exists("covar_terms")) {
  covar_terms <- "PC1 + PC2 + PC3 + PC4 + PC5 + D_MSGP1_AU + D_MSGP2_AU + D_MSPGX_AU"
}

# -------------------------- Sanity checks & prep ------------------------------
stopifnot(exists("info"), exists("dose"), exists("pheno"), exists("cov"),
          exists("allelenames"), exists("pheno_name"))

info$pos <- as.character(info$pos)
if (!"tag" %in% names(info)) info$tag <- paste0(info$pos, ":", info$AA)

# -------------------------- Helper: position predictors (FIXED) ---------------
# IMPORTANT: use different variable names (gene0, pos0) to avoid subset scoping bugs.
mk_position_predictors <- function(dat, info, gene0, pos0, allelenames) {
  sub <- info[info$gene == gene0 & info$pos == as.character(pos0),
              c("hla","AA"), drop = FALSE]
  if (nrow(sub) == 0) return(list(dat = dat, added = character(0)))
  sub <- sub[sub$hla %in% allelenames, , drop = FALSE]
  if (nrow(sub) == 0) return(list(dat = dat, added = character(0)))

  residues <- sort(unique(sub$AA))
  added <- character(0)
  for (aa in residues) {
    hlas <- sub$hla[sub$AA == aa]
    hlas_present <- intersect(hlas, colnames(dat))
    if (length(hlas_present) == 0) next
    colname <- paste0(gene0, "_p", pos0, "_", aa)
    dat[[colname]] <- if (length(hlas_present) == 1L) dat[[hlas_present]]
                      else rowSums(dat[, hlas_present, drop = FALSE], na.rm = TRUE)
    added <- c(added, colname)
  }

  # Drop the most frequent residue at that position as reference (avoid collinearity)
  if (length(added) >= 2) {
    means <- vapply(added, function(x) mean(dat[[x]], na.rm = TRUE), numeric(1))
    ref <- names(sort(means, decreasing = TRUE))[1]
    dat[[ref]] <- NULL
    added <- setdiff(added, ref)
  }
  list(dat = dat, added = added)
}

# -------------------------- Helpers: fit & Nagelkerke R² ----------------------
nagelkerke_r2 <- function(model_null, model_fit) {
  L0 <- as.numeric(logLik(model_null))
  L1 <- as.numeric(logLik(model_fit))
  N  <- nobs(model_fit)
  num <- 1 - exp((2/N) * (L0 - L1))
  den <- 1 - exp((2/N) * (L0))
  pmax(0, num / den)  # guard tiny negatives
}

fit_model_getR2 <- function(dat, covar_terms, pheno_name, add_predictors = character(0), model_null = NULL) {
  cols <- c(add_predictors)
  rhs  <- if (length(cols)) paste(c(paste(make.names(cols), collapse = " + "), covar_terms), collapse = " + ")
          else covar_terms
  f_null <- as.formula(paste0(pheno_name, " ~ ", covar_terms))
  f_fit  <- as.formula(paste0(pheno_name, " ~ ", rhs))
  dat <- `colnames<-`(dat, make.names(colnames(dat)))

  m0 <- if (is.null(model_null)) glm(f_null, data = dat, family = binomial()) else model_null
  m1 <- glm(f_fit, data = dat, family = binomial())
  list(R2 = nagelkerke_r2(m0, m1), m0 = m0, m1 = m1)
}

# -------------------------- Build analysis dataset ----------------------------
dat <- merge(dose, pheno, by = "IID")
dat <- merge(dat,  cov,   by = "IID")
dat[[pheno_name]] <- dat[[pheno_name]] - 1  # 0/1 as in your script

# Null (covariates-only) model once
f_null <- as.formula(paste0(pheno_name, " ~ ", covar_terms))
dat <- `colnames<-`(dat, make.names(colnames(dat)))
m_null <- glm(f_null, data = dat, family = binomial())

# -------------------------- Create predictors for all positions ---------------
added_by_pos <- list()
dat_aug <- dat
for (gene in names(pos_spec)) {
  for (pos in pos_spec[[gene]]) {
    res <- mk_position_predictors(dat_aug, info, gene, pos, allelenames)
    dat_aug <- res$dat
    added_by_pos[[paste0(gene, ":", pos)]] <- res$added
  }
}

# Optional quick diagnostics:
# print(sapply(added_by_pos, length))  # should differ by position, not all 22

# -------------------------- (1) Per-position R² (vs null) ---------------------
per_pos <- lapply(names(added_by_pos), function(k) {
  preds <- added_by_pos[[k]]
  if (length(preds) == 0) return(data.frame(key = k, R2 = NA_real_))
  out <- fit_model_getR2(dat_aug, covar_terms, pheno_name, add_predictors = preds, model_null = m_null)
  data.frame(key = k, R2 = out$R2)
})
per_pos_df <- do.call(rbind, per_pos)

# -------------------------- (2) Cumulative R² in chosen order -----------------
order_keys <- c("DRB1:71","DRB1:13","DRB1:86","A:95","DPB1:65")
cum_rows <- list(); current_preds <- character(0); prev_R2 <- 0
for (k in order_keys) {
  preds <- added_by_pos[[k]]
  if (is.null(preds) || length(preds) == 0) {
    cum_rows[[k]] <- data.frame(key = k, cum_R2 = NA_real_, delta_R2 = NA_real_, n_params_added = 0)
    next
  }
  # de-duplicate in case of overlap
  current_preds <- unique(c(current_preds, preds))
  out <- fit_model_getR2(dat_aug, covar_terms, pheno_name, add_predictors = current_preds, model_null = m_null)
  cum_R2 <- out$R2
  cum_rows[[k]] <- data.frame(key = k, cum_R2 = cum_R2, delta_R2 = cum_R2 - prev_R2, n_params_added = length(preds))
  prev_R2 <- cum_R2
}
cum_df <- do.call(rbind, cum_rows)

# -------------------------- (3) Partial/unique R² by drop-one -----------------
all_preds <- unlist(added_by_pos, use.names = FALSE)
all_preds <- unique(all_preds[!is.na(all_preds)])
full_fit <- fit_model_getR2(dat_aug, covar_terms, pheno_name, add_predictors = all_preds, model_null = m_null)
R2_full <- full_fit$R2

partial_rows <- lapply(names(added_by_pos), function(k) {
  preds_k <- added_by_pos[[k]]
  keep <- setdiff(all_preds, preds_k)
  if (length(preds_k) == 0 || length(keep) == 0)
    return(data.frame(key = k, partial_R2 = NA_real_))
  out <- fit_model_getR2(dat_aug, covar_terms, pheno_name, add_predictors = keep, model_null = m_null)
  data.frame(key = k, partial_R2 = R2_full - out$R2)
})
partial_df <- do.call(rbind, partial_rows)

# -------------------------- Merge & save --------------------------------------
summary_df <- Reduce(function(a,b) merge(a,b,by="key",all=TRUE),
                     list(per_pos_df, cum_df, partial_df))
summary_df$gene <- sub(":.*","", summary_df$key)
summary_df$pos  <- sub(".*:","", summary_df$key)
summary_df <- summary_df[, c("gene","pos","R2","cum_R2","delta_R2","partial_R2","n_params_added","key")]

summary_sort<-summary_df[order(match(summary_df$key, order_keys)), ]
out_summary <- "AA_positions_variance_explained_summary.tsv"
write.table(summary_sort, out_summary, sep="\t", quote=FALSE, row.names=FALSE)
cat("Wrote ", out_summary, "\n", sep="")

# Handy printout:
print(summary_df[order(match(summary_df$key, order_keys)), ], row.names = FALSE)
