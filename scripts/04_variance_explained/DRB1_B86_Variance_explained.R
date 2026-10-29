if (!nzchar(Sys.getenv("HLA_PROJECT_ROOT", unset = ""))) stop("Set HLA_PROJECT_ROOT to the input-data root before running this historical script.")
# =============================================================================
# DRB1 beta86 variance explained analysis
# 1) Primary: DRB1 beta86 position-level test conditional on HLA_DRB1_15_01
# 2) Supporting: explicit beta86V vs beta86G model
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

options(stringsAsFactors = FALSE)
`%&%` <- function(a, b) paste0(a, b)

# -------------------------- Toggle: minimal loader ----------------------------
RUN_MINIMAL_LOADER <- TRUE

if (RUN_MINIMAL_LOADER) {
  data_assoc  <- file.path(Sys.getenv("HLA_PROJECT_ROOT", unset = ""), "project/HLA_Imputation/HLA_analyses_tutorial/data_assoc/")
  raw_prefix  <- "AU_HLA_imputed_merged.keep.qc_filtered.QCed_HLA_tf.imputed"
  imputed     <- "AU_HLA_imputed_merged.keep.qc_filtered.QCed_HLA_tf"
  pheno_file  <- "Case_control.txt"
  cov_name    <- "Case_Control_covariates.txt"
  pheno_name  <- "MS_pheno"

  info <- readRDS(data_assoc %&% "HLA_DICTIONARY_AA.hg19.imgt3320.AA_tf.in_ref.rds")

  dose <- read.table(raw_prefix %&% ".raw", header = TRUE, check.names = FALSE)
  dose <- dose[, c(-1, -3:-6)]

  allelenames <- as.character(read.table(imputed %&% ".allele_name.txt")$V1)
  colnames(dose) <- c("IID", allelenames)

  pheno <- read.table(pheno_file, header = TRUE)
  cov   <- read.table(cov_name, header = TRUE)
  cov   <- cov[, colnames(cov) %in% c(
    "IID", "PC1", "PC2", "PC3", "PC4", "PC5",
    "D_MSGP1_AU", "D_MSGP2_AU", "D_MSPGX_AU"
  )]
}

# =============================================================================
# User settings
# =============================================================================

gene0 <- "DRB1"
pos0  <- "86"

condition_alleles <- c("HLA_DRB1_15_01")

covar_terms <- "PC1 + PC2 + PC3 + PC4 + PC5 + D_MSGP1_AU + D_MSGP2_AU + D_MSPGX_AU"

out_prefix <- "DRB1_B86_condition_DRB1_1501"

# =============================================================================
# Prepare data
# =============================================================================

info$pos <- as.character(info$pos)

dat <- merge(dose, pheno, by = "IID")
dat <- merge(dat, cov, by = "IID")

if (all(na.omit(unique(dat[[pheno_name]])) %in% c(1, 2))) {
  dat[[pheno_name]] <- dat[[pheno_name]] - 1
}

colnames(dat) <- make.names(colnames(dat))
allelenames_safe <- make.names(allelenames)
condition_alleles_safe <- make.names(condition_alleles)
pheno_name_safe <- make.names(pheno_name)

# =============================================================================
# Helper functions
# =============================================================================

nagelkerke_r2 <- function(model_null, model_fit) {
  L0 <- as.numeric(logLik(model_null))
  L1 <- as.numeric(logLik(model_fit))
  N  <- nobs(model_fit)

  num <- 1 - exp((2 / N) * (L0 - L1))
  den <- 1 - exp((2 / N) * L0)

  pmax(0, num / den)
}

make_formula <- function(pheno, terms) {
  terms <- terms[!is.na(terms) & terms != ""]
  as.formula(paste(pheno, "~", paste(terms, collapse = " + ")))
}

make_AA_residue_predictors <- function(dat, info, gene0, pos0, allelenames) {

  sub <- info[
    info$gene == gene0 & info$pos == as.character(pos0),
    c("hla", "AA"),
    drop = FALSE
  ]

  sub <- sub[sub$hla %in% allelenames, , drop = FALSE]

  if (nrow(sub) == 0) {
    stop("No AA dictionary rows found for ", gene0, " position ", pos0)
  }

  residues <- sort(unique(sub$AA))
  added <- character(0)

  for (aa in residues) {
    hlas <- sub$hla[sub$AA == aa]
    hlas_safe <- make.names(hlas)
    hlas_present <- intersect(hlas_safe, colnames(dat))

    if (length(hlas_present) == 0) next

    pred_name <- make.names(paste0(gene0, "_p", pos0, "_", aa))

    dat[[pred_name]] <- if (length(hlas_present) == 1) {
      dat[[hlas_present]]
    } else {
      rowSums(dat[, hlas_present, drop = FALSE], na.rm = TRUE)
    }

    added <- c(added, pred_name)
  }

  list(
    dat = dat,
    predictors = added
  )
}

drop_reference_residue <- function(dat, predictors, gene0, pos0, reference_residue = NULL) {

  if (length(predictors) < 2) {
    return(list(
      dat = dat,
      kept_predictors = predictors,
      dropped_predictor = NA_character_,
      dropped_residue = NA_character_
    ))
  }

  if (is.null(reference_residue)) {
    means <- vapply(predictors, function(x) mean(dat[[x]], na.rm = TRUE), numeric(1))
    dropped_predictor <- names(sort(means, decreasing = TRUE))[1]
  } else {
    dropped_predictor <- make.names(paste0(gene0, "_p", pos0, "_", reference_residue))

    if (!(dropped_predictor %in% predictors)) {
      stop("Requested reference residue not found: ", reference_residue)
    }
  }

  dropped_residue <- sub(make.names(paste0(gene0, "_p", pos0, "_")), "", dropped_predictor)

  dat[[dropped_predictor]] <- NULL
  kept_predictors <- setdiff(predictors, dropped_predictor)

  list(
    dat = dat,
    kept_predictors = kept_predictors,
    dropped_predictor = dropped_predictor,
    dropped_residue = dropped_residue
  )
}

# =============================================================================
# Create beta86 residue predictors
# =============================================================================

aa_obj <- make_AA_residue_predictors(
  dat = dat,
  info = info,
  gene0 = gene0,
  pos0 = pos0,
  allelenames = allelenames
)

dat_aa <- aa_obj$dat
all_b86_predictors <- aa_obj$predictors

cat("\nCreated beta86 residue predictors:\n")
print(all_b86_predictors)

# =============================================================================
# 1) PRIMARY RESULT:
# DRB1 beta86 position-level test conditional on HLA_DRB1_15_01
# =============================================================================

# For position-level test, drop the most frequent residue as reference.
# This tests the whole beta86 position using K-1 degrees of freedom.

primary_ref <- drop_reference_residue(
  dat = dat_aa,
  predictors = all_b86_predictors,
  gene0 = gene0,
  pos0 = pos0,
  reference_residue = NULL
)

dat_primary <- primary_ref$dat
b86_position_terms <- primary_ref$kept_predictors

base_terms <- c(covar_terms, condition_alleles_safe)

primary_reduced_formula <- make_formula(
  pheno = pheno_name_safe,
  terms = base_terms
)

primary_full_formula <- make_formula(
  pheno = pheno_name_safe,
  terms = c(base_terms, b86_position_terms)
)

primary_reduced_model <- glm(primary_reduced_formula, data = dat_primary, family = binomial())
primary_full_model    <- glm(primary_full_formula,    data = dat_primary, family = binomial())

primary_R2 <- nagelkerke_r2(primary_reduced_model, primary_full_model)
primary_lrt <- anova(primary_reduced_model, primary_full_model, test = "Chisq")

primary_summary <- data.frame(
  analysis = "Primary_position_level_test",
  gene = gene0,
  position = pos0,
  conditioned_on = paste(condition_alleles, collapse = ","),
  dropped_reference_residue = primary_ref$dropped_residue,
  tested_terms = paste(b86_position_terms, collapse = ","),
  reduced_formula = paste(deparse(primary_reduced_formula), collapse = ""),
  full_formula = paste(deparse(primary_full_formula), collapse = ""),
  Nagelkerke_R2_increment = primary_R2,
  LRT_deviance = primary_lrt$Deviance[2],
  LRT_df = primary_lrt$Df[2],
  LRT_p = primary_lrt$`Pr(>Chi)`[2],
  stringsAsFactors = FALSE
)

write.table(
  primary_summary,
  paste0(out_prefix, "_PRIMARY_B86_position_level.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n===== PRIMARY RESULT: DRB1 beta86 position-level test =====\n")
print(primary_summary, row.names = FALSE)

# =============================================================================
# 2) SUPPORTING RESULT:
# Explicit beta86V vs beta86G model
# =============================================================================

# To explicitly model beta86V, drop G as the reference.
# Model becomes:
# MS_pheno ~ covariates + HLA_DRB1_15_01 + DRB1_p86_V

support_ref <- drop_reference_residue(
  dat = dat_aa,
  predictors = all_b86_predictors,
  gene0 = gene0,
  pos0 = pos0,
  reference_residue = "G"
)

dat_support <- support_ref$dat
support_terms <- support_ref$kept_predictors

if (!any(grepl("_V$", support_terms))) {
  warning("DRB1_p86_V was not retained. Check available beta86 residues.")
}

support_reduced_formula <- make_formula(
  pheno = pheno_name_safe,
  terms = base_terms
)

support_full_formula <- make_formula(
  pheno = pheno_name_safe,
  terms = c(base_terms, support_terms)
)

support_reduced_model <- glm(support_reduced_formula, data = dat_support, family = binomial())
support_full_model    <- glm(support_full_formula,    data = dat_support, family = binomial())

support_R2 <- nagelkerke_r2(support_reduced_model, support_full_model)
support_lrt <- anova(support_reduced_model, support_full_model, test = "Chisq")

support_coef <- summary(support_full_model)$coefficients

support_coef_df <- as.data.frame(support_coef)
support_coef_df$term <- rownames(support_coef_df)
rownames(support_coef_df) <- NULL

support_coef_df <- support_coef_df[
  support_coef_df$term %in% support_terms,
]

support_coef_df$OR <- exp(support_coef_df$Estimate)
support_coef_df$CI_lower <- exp(support_coef_df$Estimate - 1.96 * support_coef_df$`Std. Error`)
support_coef_df$CI_upper <- exp(support_coef_df$Estimate + 1.96 * support_coef_df$`Std. Error`)

support_summary <- data.frame(
  analysis = "Supporting_explicit_beta86V_vs_beta86G",
  gene = gene0,
  position = pos0,
  conditioned_on = paste(condition_alleles, collapse = ","),
  dropped_reference_residue = support_ref$dropped_residue,
  tested_terms = paste(support_terms, collapse = ","),
  reduced_formula = paste(deparse(support_reduced_formula), collapse = ""),
  full_formula = paste(deparse(support_full_formula), collapse = ""),
  Nagelkerke_R2_increment = support_R2,
  LRT_deviance = support_lrt$Deviance[2],
  LRT_df = support_lrt$Df[2],
  LRT_p = support_lrt$`Pr(>Chi)`[2],
  stringsAsFactors = FALSE
)

write.table(
  support_summary,
  paste0(out_prefix, "_SUPPORTING_B86V_vs_B86G_model.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  support_coef_df,
  paste0(out_prefix, "_SUPPORTING_B86V_vs_B86G_coefficients.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n===== SUPPORTING RESULT: explicit beta86V vs beta86G =====\n")
print(support_summary, row.names = FALSE)

cat("\n===== SUPPORTING COEFFICIENTS =====\n")
print(support_coef_df, row.names = FALSE)

# =============================================================================
# Interpretation helper
# =============================================================================

cat("\nInterpretation:\n")
cat("Primary test:\n")
cat("  Tests whether DRB1 beta86 explains additional variance after conditioning on HLA_DRB1_15_01.\n")
cat("  This is the position-level test and should be the main result.\n\n")

cat("Supporting test:\n")
cat("  Uses G as reference so that DRB1_p86_V is explicitly included in the model.\n")
cat("  The beta/OR for DRB1_p86_V represents the effect of beta86V relative to beta86G,\n")
cat("  conditional on HLA_DRB1_15_01 and covariates.\n")
