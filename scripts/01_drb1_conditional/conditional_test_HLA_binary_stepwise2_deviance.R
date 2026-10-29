# ------------------------------- Libraries ------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(tidyr)
})

`%&%` <- function(a,b) paste0(a,b)

args <- commandArgs(TRUE)

## ----- CONFIG -----
imputed     <- args[1]
raw_prefix  <- args[2]
pheno_file  <- args[3]   # e.g., MS_median_scores_keep.txt
pheno_name  <- args[4]
cov_name    <- args[5]
HLA         <- args[6]   # e.g., "DRB1"
# NEW: optional 7th arg for positions to ignore (e.g., "96" or "DRB1:96, DRB1:13")
ignore_arg  <- if (length(args) >= 7) args[7] else ""

# In round 1, select the most significant amino-acid position, consistent with
# the omnibus analysis. In conditional rounds, select the candidate producing
# the largest change in deviance. Retain a selected position only if its
# likelihood-ratio test P value is below the study-wide threshold.
covar_terms_default <- "PC1 + PC2 + PC3 + PC4 + PC5 + D_MSGP1_AU + D_MSGP2_AU + D_MSPGX_AU"

covar_terms_round3 <- covar_terms_default

alpha <- 1e-5

# ------------------------------- Paths ----------------------------------------
data_assoc <- Sys.getenv("HLA_DICTIONARY_DIR", unset = "")
if (!nzchar(data_assoc)) stop("Set HLA_DICTIONARY_DIR to the directory containing the HLA dictionary")
data_assoc <- paste0(normalizePath(data_assoc, mustWork = TRUE), "/")

# ------------------------------ Load dictionary -------------------------------
info <- readRDS(data_assoc %&% "HLA_DICTIONARY_AA.hg19.imgt3320.AA_tf.in_ref.rds")
#info <- read.table(data_assoc %&% "HLA_DICTIONARY_AA.hg19.imgt3320.AA_tf.in_ref_original.txt",header = TRUE,stringsAsFactors = FALSE)
info <- info[info$gene == HLA, , drop = FALSE]
info$pos <- as.character(info$pos)  # keep positions as characters
info$pos_lab <- paste0("p", info$pos)
info$tag <- paste0(info$pos, ":", info$AA)

# ---------------------------- Apply ignore filter -----------------------------
ignore_tokens <- if (nzchar(ignore_arg)) unlist(strsplit(ignore_arg, "[,;\\s]+")) else character(0)
ignore_pos <- unique(gsub("^.*?:", "", ignore_tokens))   # strip any "DRB1:" prefixes
ignore_pos <- ignore_pos[nzchar(ignore_pos)]

if (length(ignore_pos) > 0) {
  info <- info[!(info$pos %in% ignore_pos), , drop = FALSE]
  message("Ignoring positions for ", HLA, ": ", paste(ignore_pos, collapse = ", "))
}

# --------------------------- Build position universe --------------------------
allpos <- sort(unique(info$pos))
if (length(allpos) == 0) stop("No positions left to test after applying ignore list.")

# --------------------- wide allele → AA table (by position) -------------------
allele_res_wide <- info %>%
  select(hla, pos, pos_lab, AA) %>%
  distinct() %>%
  pivot_wider(names_from = pos_lab, values_from = AA)

# -------------------------- Load dosages / names ------------------------------
dose <- read.table(raw_prefix %&% ".raw", header = TRUE, check.names = FALSE)
dose <- dose[, c(-1, -3:-6)]  # drop FID, PAT, MAT, SEX, PHENOTYPE
allelenames <- as.character(read.table(imputed %&% ".QCed_HLA_tf.allele_name.txt")$V1)
colnames(dose) <- c("IID", allelenames)

# keep only alleles with dosages in the wide table
allele_res_wide <- allele_res_wide %>% filter(hla %in% allelenames)

# -------------------------- phenotype + covariates -----------------------------
pheno <- read.table(pheno_file, header = TRUE)
cov   <- read.table(cov_name, header = TRUE)
cov   <- cov[, (colnames(cov) %in% c("IID","PC1","PC2","PC3","PC4","PC5","D_MSGP1_AU","D_MSGP2_AU","D_MSPGX_AU"))]

dat <- merge(dose, pheno, by = "IID")
dat <- merge(dat, cov, by = "IID")
dat[[pheno_name]] <- dat[[pheno_name]] - 1  # make binary 0/1

# ------------------------------- Helpers --------------------------------------
# Ensure tag column present
ensure_info_tags <- function(info) {
  info$pos <- as.character(info$pos)
  if (!"tag" %in% names(info)) info$tag <- paste0(info$pos, ":", info$AA)
  info
}
info <- ensure_info_tags(info)

# Pre-split HLA lists by tag (AFTER ignore filter!)
hla_by_tag <- split(info$hla, info$tag)

# ---- Builder (wide table): combine residues at locked/candidate into hap keys
make_group_map <- function(allele_res_wide, locked = character(0), candidate_pos = NULL) {
  positions_chr <- unique(as.character(c(locked, candidate_pos)))
  if (length(positions_chr) == 0L) return(NULL)

  suppressWarnings(num <- as.numeric(positions_chr))
  ord <- if (all(!is.na(num))) order(num) else order(positions_chr)
  positions_chr <- positions_chr[ord]

  cols <- paste0("p", positions_chr)
  missing <- setdiff(cols, colnames(allele_res_wide))
  if (length(missing) > 0L) return(NULL)

  tab <- as.data.frame(allele_res_wide[, c("hla", cols), drop = FALSE])
  if (nrow(tab) == 0L) return(NULL)

  keep <- stats::complete.cases(tab[, cols, drop = FALSE])
  if (!any(keep)) return(NULL)
  tab <- tab[keep, , drop = FALSE]

  key <- apply(tab[, cols, drop = FALSE], 1, function(v) {
    paste(paste0(positions_chr, ":", v), collapse = "_")
  })

  data.frame(
    hap = key,
    hla = tab$hla,
    pos = if (length(positions_chr)) tail(positions_chr, 1) else NA_character_,
    stringsAsFactors = FALSE
  )
}

# ---- Builder (info/tag-based): exact "xtag_ytag" style via intersect on tags
make_res_map <- function(info, hla_by_tag, locked = character(0), candidate_pos = NULL) {
  positions <- unique(as.character(c(locked, candidate_pos)))
  if (length(positions) == 0L) return(NULL)

  suppressWarnings(num <- as.numeric(positions))
  ord <- if (all(!is.na(num))) order(num) else order(positions)
  positions <- positions[ord]

  tags_list <- lapply(positions, function(p) unique(info$tag[info$pos == p]))
  if (any(vapply(tags_list, length, integer(1)) == 0L)) return(NULL)

  comb <- do.call(expand.grid, c(tags_list, stringsAsFactors = FALSE))
  names(comb) <- paste0("tag_", positions)

  out <- vector("list", nrow(comb))
  j <- 0L
  for (i in seq_len(nrow(comb))) {
    tags_vec <- as.character(unlist(comb[i, ], use.names = FALSE))
    hla_sets <- lapply(tags_vec, function(tg) hla_by_tag[[tg]])
    if (any(vapply(hla_sets, is.null, logical(1)))) next
    inter <- Reduce(intersect, hla_sets)
    if (length(inter) == 0L) next

    hap_key <- paste(tags_vec, collapse = "_")
    j <- j + 1L
    out[[j]] <- data.frame(
      hap = hap_key,
      hla = inter,
      pos = if (length(positions)) tail(positions, 1) else NA_character_,
      stringsAsFactors = FALSE
    )
  }
  if (j == 0L) return(NULL)
  res <- do.call(rbind, out[seq_len(j)])
  rownames(res) <- NULL
  res
}

# ---- Add grouped dosages to dat for each hap group
add_group_dosages <- function(dat, group_map, allelenames) {
  adopted <- character(0)
  if (is.null(group_map) || nrow(group_map) == 0L) return(list(adopted = adopted, dat = dat))

  gm <- group_map[group_map$hla %in% allelenames, , drop = FALSE]
  if (nrow(gm) == 0L) return(list(adopted = adopted, dat = dat))

  for (h in unique(gm$hap)) {
    hlas <- gm$hla[gm$hap == h]
    hlas_present <- intersect(hlas, colnames(dat))
    if (length(hlas_present) == 0L) next

    dat[[h]] <- if (length(hlas_present) == 1L) {
      dat[[hlas_present]]
    } else {
      rowSums(dat[, hlas_present, drop = FALSE], na.rm = TRUE)
    }
    adopted <- c(adopted, h)
  }
  list(adopted = adopted, dat = dat)
}

# ---- Compare reduced vs full models (logistic)
compare_models <- function(dat, adopted_reduced, adopted_full, covar_terms) {
  rhs_reduced <- if (length(adopted_reduced) > 0) {
    paste(c(paste(make.names(adopted_reduced), collapse = " + "), covar_terms), collapse = " + ")
  } else {
    covar_terms
  }
  rhs_full <- paste(c(paste(make.names(adopted_full), collapse = " + "), covar_terms), collapse = " + ")

  f_red  <- as.formula(paste0(pheno_name, " ~ ", rhs_reduced))
  f_full <- as.formula(paste0(pheno_name, " ~ ", rhs_full))
  print(f_full)
  colnames(dat) <- make.names(colnames(dat))
  m_red  <- glm(f_red,  data = dat, family = binomial(link="logit"))
  m_full <- glm(f_full, data = dat, family = binomial(link="logit"))
  a <- anova(m_red, m_full, test = "Chisq")

  # For nested models, the likelihood-ratio-test degrees of freedom equal the
  # difference in the ranks of the two fitted model matrices. This can be less
  # than the nominal difference in the number of haplotype groups when columns
  # are absent or linearly dependent.
  rank_reduced <- m_red$rank
  rank_full    <- m_full$rank
  rank_diff    <- rank_full - rank_reduced

  if (!isTRUE(all.equal(as.numeric(a$Df[2]), as.numeric(rank_diff)))) {
    warning(
      "LRT degrees of freedom (", a$Df[2],
      ") do not match the model-rank difference (", rank_diff, ")."
    )
  }

  list(
    deviance  = a$Deviance[2],
    Df        = a$Df[2],
    pval      = a$`Pr(>Chi)`[2],
    Res.DF    = a$`Resid. Df`[2],
    Resid_Dev = a$`Resid. Dev`[2],
    Rank_Reduced = rank_reduced,
    Rank_Full    = rank_full,
    Rank_Diff    = rank_diff
  )
}

# ------------------------------- Stepwise loop --------------------------------
remaining <- allpos
locked    <- character(0)
round_summaries <- list()
cond_pos <- NA

repeat {
  adopted_red <- character(0)

  # Use the special covariate set only in round 3, when the model is
  # conditional on amino-acid positions 71 and 13. Other rounds retain the
  # default covariate specification.
  is_round3_71_13 <-
    length(locked) == 2L &&
    setequal(as.character(locked), c("71", "13"))

  covar_terms_active <- if (is_round3_71_13) {
    covar_terms_round3
  } else {
    covar_terms_default
  }
  
  if (length(locked) > 0) {
    cond_pos <- paste(locked, collapse = ":")
    gm_red   <- make_res_map(info, hla_by_tag, locked = locked)
    tmp_red  <- add_group_dosages(dat, gm_red, allelenames)
    dat      <- tmp_red$dat
    adopted_red <- tmp_red$adopted
  }

  if (length(remaining) == 0) break

  test_rows <- vector("list", length(remaining))
  irow <- 0L

  for (pos in remaining) {
    gm_full      <- make_res_map(info, hla_by_tag, locked = locked, candidate_pos = pos)
    tmp_full     <- add_group_dosages(dat, gm_full, allelenames)
    dat_full     <- tmp_full$dat
    adopted_full <- tmp_full$adopted

    irow <- irow + 1L
    if (length(adopted_full) == 0) {
      test_rows[[irow]] <- data.frame(
        HLA = HLA,
        round = length(locked) + 1,
        cond_pos = if (length(locked)) cond_pos else NA,
        add_pos  = pos,
        deviance = NA_real_,
        Df       = NA_real_,
        pval     = NA_real_,
        Res.DF   = NA_real_,
        Resid_Dev= NA_real_,
        Rank_Reduced = NA_real_,
        Rank_Full    = NA_real_,
        Rank_Diff    = NA_real_,
        Covariates   = covar_terms_active,
        stringsAsFactors = FALSE
      )
      next
    }
    cmp <- compare_models(
      dat = dat_full,
      adopted_reduced = adopted_red,
      adopted_full = adopted_full,
      covar_terms = covar_terms_active
    )
    test_rows[[irow]] <- data.frame(
      HLA = HLA,
      round = length(locked) + 1,
      cond_pos = if (length(locked)) cond_pos else NA,
      add_pos  = pos,
      deviance = cmp$deviance,
      Df       = cmp$Df,
      pval     = cmp$pval,
      Res.DF   = cmp$Res.DF,
      Resid_Dev= cmp$Resid_Dev,
      Rank_Reduced = cmp$Rank_Reduced,
      Rank_Full    = cmp$Rank_Full,
      Rank_Diff    = cmp$Rank_Diff,
      Covariates   = covar_terms_active,
      stringsAsFactors = FALSE
    )
  }

  # Round 1 follows the omnibus analysis and selects the smallest P value.
  # Conditional rounds select the largest additional improvement in deviance.
  # In every round, the selected position is retained only if P < alpha.
  round_df <- bind_rows(test_rows)
  if (length(locked) == 0L) {
    round_df <- round_df %>% arrange(pval, desc(deviance))
  } else {
    round_df <- round_df %>% arrange(desc(deviance), pval)
  }
  round_summaries[[length(locked) + 1]] <- round_df

  if (nrow(round_df) == 0) break
  best <- round_df[1, ]
  if (is.na(best$pval) || best$pval >= alpha) break

  locked    <- c(locked, best$add_pos)
  remaining <- setdiff(remaining, best$add_pos)
}

# --------------------------------- Results ------------------------------------
summary_all_rounds <- bind_rows(round_summaries)
locked_positions_in_order <- locked

cat("\n==== Stepwise conditional haplotype test complete ====\n")
cat("HLA gene: ", HLA, "\n", sep = "")
cat("Selection rule: minimum P value in round 1; maximum change in deviance in conditional rounds\n")
cat("Retention threshold: conditional LRT P < ", format(alpha, scientific = TRUE), "\n", sep = "")
cat("Selected positions (in order): ", if (length(locked_positions_in_order)) paste(locked_positions_in_order, collapse = " -> ") else "none", "\n", sep = "")
if (length(ignore_pos) > 0) cat("Ignored positions: ", paste(ignore_pos, collapse = ", "), "\n", sep = "")
cat("\nTop of results table:\n")
print(utils::head(summary_all_rounds, 10))

out_file <- "conditional_test_HLA_" %&% HLA %&% "_" %&% pheno_name %&% "_all.txt"
write.table(summary_all_rounds, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nWrote: ", out_file, "\n", sep = "")
