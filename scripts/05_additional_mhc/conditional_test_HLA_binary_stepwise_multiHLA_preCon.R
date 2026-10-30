#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(tidyr)
})

`%&%` <- function(a,b) paste(a,b, sep = "")

args <- commandArgs(TRUE)

## ===================== CONFIG =====================
imputed     <- args[1]  # e.g. AU_HLA_imputed_merged.keep.qc_filtered
raw_prefix  <- args[2]  # same prefix (has *.raw)
pheno_file  <- args[3]  # e.g. Case_control.txt
pheno_name  <- args[4]  # e.g. MS_pheno
cov_name    <- args[5]  # e.g. Case_Control_covariates.txt
HLA_arg     <- args[6]  # e.g. "DRB1,DQA1,DQB1"

## NEW(7th): prelocked positions via CLI. Example: "DRB1:71,DRB1:13,DRB1:37,DRB1:86"
prelocked_arg <- if (length(args) >= 7) args[7] else ""
parse_list <- function(s) {
  if (is.null(s) || nchar(trimws(s)) == 0 || s %in% c("-", "NA", "none", "None")) return(character(0))
  out <- unlist(strsplit(s, "[,;\\s]+"))
  trimws(out[nzchar(out)])
}
PRELOCKED_POS <- parse_list(prelocked_arg)

## NEW(8th): positions to ignore. Example: "B:326" or "B:326;A:62"
ignore_arg <- if (length(args) >= 8) args[8] else ""
IGNORE_POS <- parse_list(ignore_arg)
print(IGNORE_POS)

covar_terms <- "PC1 + PC2 + PC3 + PC4 + PC5"
alpha       <- 1e-5
max_rounds  <- 50

# ===============================================
# ====== REGION FILTER (now OFF by default) ======
GROOVE_ONLY <- FALSE     # set TRUE to restrict to groove
GROOVE_MIN  <- 1
GROOVE_MAX  <- 94

# Paths
data_assoc <- Sys.getenv("HLA_DICTIONARY_DIR", unset = "")
if (!nzchar(data_assoc)) stop("Set HLA_DICTIONARY_DIR to the directory containing the HLA dictionary")
data_assoc <- paste0(sub("/+$", "", data_assoc), "/")
## ==================================================

# ----- parse genes -----
HLA_GENES <- strsplit(HLA_arg, ",")[[1]] %>% trimws()
HLA_LABEL <- paste(HLA_GENES, collapse = ",")

# ----- load dictionary -----
info_all <- readRDS(data_assoc %&% "HLA_DICTIONARY_AA.hg19.imgt3320.AA_tf.in_ref.rds")

info <- info_all %>%
  filter(gene %in% HLA_GENES) %>%
  mutate(
    pos_chr = as.character(pos),
    pos_num = suppressWarnings(as.numeric(pos_chr))
  )
if (nrow(info) == 0) stop("No AA dictionary rows for requested genes: ", HLA_LABEL)

# (Optional) groove filter for *candidates*; prelocked will be added back explicitly.
if (GROOVE_ONLY) {
  classII <- c("DRB1","DQA1","DQB1")
  info <- info %>%
    filter(
      (!gene %in% classII) |
      (gene %in% classII & !is.na(pos_num) & pos_num >= GROOVE_MIN & pos_num <= GROOVE_MAX)
    )
}

# Ensure PRELOCKED rows are present even if groove filter would exclude them
if (length(PRELOCKED_POS) > 0) {
  prelocked_rows <- info_all %>%
    mutate(pos_chr = as.character(pos)) %>%
    filter(paste0(gene, ":", pos_chr) %in% PRELOCKED_POS) %>%
    mutate(pos_num = suppressWarnings(as.numeric(pos_chr)))
  if (nrow(prelocked_rows) > 0) {
    info <- bind_rows(
      info,
      anti_join(prelocked_rows,
                info %>% mutate(pos_chr = as.character(pos)),
                by = c("gene","hla","AA","pos"))
    )
  }
}

# Build IDs BEFORE ignore filtering (so we can resolve overlaps cleanly)
info <- info %>%
  mutate(
    pos       = as.character(pos),
    pos_id    = paste0(gene, ":", pos),         # e.g. DRB1:71 or B:326
    tag       = paste0(gene, ":", pos, ":", AA) # e.g. DRB1:71:R
  )

# If any position is both prelocked and ignored, keep it (prelocked wins) and warn
if (length(PRELOCKED_POS) > 0 && length(IGNORE_POS) > 0) {
  overlap <- intersect(PRELOCKED_POS, IGNORE_POS)
  if (length(overlap) > 0) {
    warning("Positions appear in both PRELOCKED and IGNORE lists; keeping them locked: ",
            paste(overlap, collapse = ", "))
    IGNORE_POS <- setdiff(IGNORE_POS, overlap)
  }
}

# ---------------------------- Apply ignore filter -----------------------------
if (length(IGNORE_POS) > 0) {
  info <- info %>% filter(!(pos_id %in% IGNORE_POS))
  message("Ignoring positions: ", paste(IGNORE_POS, collapse = ", "))
}

# tag -> alleles (AFTER ignore)
hla_by_tag <- split(info$hla, info$tag)

ALL_POS <- sort(unique(info$pos_id))
if (length(ALL_POS) == 0) stop("No testable positions after filtering.")

# Sanity: keep only prelocked that actually exist in dictionary after filtering
if (length(PRELOCKED_POS) > 0) {
  PRELOCKED_POS <- intersect(PRELOCKED_POS, ALL_POS)
  if (length(PRELOCKED_POS) == 0) {
    warning("None of the requested pre-locked positions found post-filtering; proceeding without prelocks.")
  }
}

# ----- data: dosages, pheno, covars -----
dose <- read.table(raw_prefix %&% ".raw", header = TRUE, check.names = FALSE)
dose <- dose[, c(-1, -3:-6)] # keep IID + alleles
allelenames <- as.character(read.table(imputed %&% ".QCed_HLA_tf.allele_name.txt")$V1)
colnames(dose) <- c("IID", allelenames)

pheno <- read.table(pheno_file, header = TRUE, check.names = FALSE)
cov   <- read.table(cov_name,   header = TRUE, check.names = FALSE)
cov   <- cov[, (colnames(cov) %in% c("IID","PC1","PC2","PC3","PC4","PC5")), drop = FALSE]

dat <- merge(dose, pheno, by = "IID")
dat <- merge(dat,  cov,   by = "IID")

# Binary phenotype 0/1 safeguard
if (is.numeric(dat[[pheno_name]]) &&
    max(dat[[pheno_name]], na.rm = TRUE) == 2 &&
    min(dat[[pheno_name]], na.rm = TRUE) == 1) {
  dat[[pheno_name]] <- dat[[pheno_name]] - 1
}
if (!all(dat[[pheno_name]] %in% c(0,1))) {
  stop("Phenotype column '", pheno_name, "' must be 0/1 or 1/2.")
}

# ----- helpers (additive per-position encoding) -----
tags_for_position <- function(info, pos_id) {
  unique(info$tag[info$pos_id == pos_id])
}

add_residue_columns <- function(dat, tags, allelenames) {
  added <- character(0)
  for (tg in tags) {
    if (tg %in% colnames(dat)) { added <- c(added, tg); next }
    hlas <- hla_by_tag[[tg]]
    if (is.null(hlas)) next
    present <- intersect(hlas, colnames(dat))
    if (length(present) == 0L) next
    dat[[tg]] <- if (length(present) == 1L) dat[[present]] else rowSums(dat[, present, drop = FALSE], na.rm = TRUE)
    added <- c(added, tg)
  }
  list(dat = dat, added = added)
}

choose_ref_tag <- function(dat, pos_tags, pheno_name) {
  avail <- intersect(pos_tags, colnames(dat))
  if (length(avail) == 0) return(NA_character_)
  idx_ctrl <- which(dat[[pheno_name]] == 0)
  means <- if (length(idx_ctrl) > 0) {
    vapply(avail, function(cn) mean(dat[idx_ctrl, cn], na.rm = TRUE), numeric(1))
  } else {
    vapply(avail, function(cn) mean(dat[[cn]], na.rm = TRUE), numeric(1))
  }
  names(which.max(means))
}

ensure_union_terms <- function(dat, info, pos_ids, allelenames, pheno_name) {
  term_names <- character(0)
  ref_map    <- setNames(rep(NA_character_, length(pos_ids)), pos_ids)
  for (pid in pos_ids) {
    pos_tags <- tags_for_position(info, pid)
    if (length(pos_tags) == 0) next
    tmp <- add_residue_columns(dat, pos_tags, allelenames)
    dat <- tmp$dat
    ref_tag <- choose_ref_tag(dat, pos_tags, pheno_name)
    ref_map[pid] <- ref_tag
    use_tags <- setdiff(intersect(pos_tags, colnames(dat)), ref_tag)
    term_names <- c(term_names, use_tags)
  }
  list(dat = dat, terms = unique(term_names), ref_map = ref_map)
}

ensure_candidate_terms <- function(dat, info, pos_id, allelenames, pheno_name) {
  pos_tags <- tags_for_position(info, pos_id)
  if (length(pos_tags) == 0) return(list(dat = dat, terms = character(0), ref_tag = NA_character_))
  tmp <- add_residue_columns(dat, pos_tags, allelenames)
  dat <- tmp$dat
  ref_tag <- choose_ref_tag(dat, pos_tags, pheno_name)
  use_tags <- setdiff(intersect(pos_tags, colnames(dat)), ref_tag)
  list(dat = dat, terms = unique(use_tags), ref_tag = ref_tag)
}

lrt_compare <- function(dat, red_terms, full_terms, covar_terms) {
  rhs_reduced <- if (length(red_terms) > 0) {
    paste(c(paste(make.names(red_terms), collapse = " + "), covar_terms), collapse = " + ")
  } else covar_terms
  rhs_full <- paste(c(paste(make.names(full_terms), collapse = " + "), covar_terms), collapse = " + ")
  f_red  <- as.formula(paste0(pheno_name, " ~ ", rhs_reduced))
  f_full <- as.formula(paste0(pheno_name, " ~ ", rhs_full))

  colnames(dat) <- make.names(colnames(dat))
  m_red  <- glm(f_red,  data = dat, family = binomial(link = "logit"))
  m_full <- glm(f_full, data = dat, family = binomial(link = "logit"))
  a <- anova(m_red, m_full, test = "Chisq")
  list(
    deviance  = a$Deviance[2],
    Df        = a$Df[2],
    pval      = a$`Pr(>Chi)`[2],
    Res.DF    = a$`Resid. Df`[2],
    Resid_Dev = a$`Resid. Dev`[2]
  )
}

# ---------------- stepwise (additive per-position) ----------------

# Pre-lock from CLI (can be empty)
locked <- PRELOCKED_POS

# Test only positions OUTSIDE DRB1
remaining <- setdiff(ALL_POS[!grepl("^DRB1:", ALL_POS)], PRELOCKED_POS)

# Also ensure ignored positions are not in remaining (redundant, but explicit)
if (length(IGNORE_POS) > 0) {
  remaining <- setdiff(remaining, IGNORE_POS)
}

round_summaries <- list()
round_index <- 0L

while (length(remaining) > 0 && length(locked) - length(PRELOCKED_POS) < max_rounds) {
  round_index <- round_index + 1L

  red <- ensure_union_terms(dat, info, locked, allelenames, pheno_name)
  dat <- red$dat
  red_terms <- red$terms

  test_rows <- vector("list", length(remaining))
  names(test_rows) <- remaining

  for (pid in remaining) {
    cand <- ensure_candidate_terms(dat, info, pid, allelenames, pheno_name)
    dat_full   <- cand$dat
    cand_terms <- cand$terms

    if (length(cand_terms) == 0) {
      test_rows[[pid]] <- data.frame(
        HLA = HLA_LABEL, round = round_index,
        cond_pos = if (length(locked)) paste(locked, collapse = ",") else NA,
        add_pos  = pid,
        deviance = NA_real_, Df = NA_real_, pval = NA_real_,
        Res.DF = NA_real_, Resid_Dev = NA_real_
      )
      next
    }

    full_terms <- unique(c(red_terms, cand_terms))
    cmp <- lrt_compare(dat_full, red_terms, full_terms, covar_terms)
    test_rows[[pid]] <- data.frame(
      HLA = HLA_LABEL, round = round_index,
      cond_pos = if (length(locked)) paste(locked, collapse = ",") else NA,
      add_pos  = pid,
      deviance = cmp$deviance, Df = cmp$Df, pval = cmp$pval,
      Res.DF = cmp$Res.DF, Resid_Dev = cmp$Resid_Dev
    )
  }

  round_df <- bind_rows(test_rows) %>% arrange(pval)
  round_summaries[[length(round_summaries) + 1L]] <- round_df

  if (nrow(round_df) == 0) break
  best <- round_df[1, ]
  if (is.na(best$pval) || best$pval >= alpha) break

  locked    <- c(locked, best$add_pos)
  remaining <- setdiff(remaining, best$add_pos)
}

# ----- results -----
summary_all_rounds <- bind_rows(round_summaries)
cat("\n==== Conditional AA test outside DRB1 complete ====\n")
cat("Genes tested: ", HLA_LABEL, "\n", sep = "")
cat("Groove-only for candidates: ", GROOVE_ONLY,
    if (GROOVE_ONLY) paste0(" (", GROOVE_MIN, "-", GROOVE_MAX, ")") else "", "\n", sep = "")
cat("Prelocked positions (from args): ",
    if (length(PRELOCKED_POS)) paste(PRELOCKED_POS, collapse = " -> ") else "<none>", "\n", sep = "")
cat("Ignored positions: ",
    if (length(IGNORE_POS)) paste(IGNORE_POS, collapse = ", ") else "<none>", "\n", sep = "")
cat("Selected (including prelocked): ", paste(locked, collapse = " -> "), "\n\n", sep = "")
suppressWarnings(print(utils::head(summary_all_rounds)))

out_fn <- "conditional_test_outsideDRB1_" %&% gsub(",", "_", HLA_LABEL) %&% "_" %&% pheno_name %&% "_additive.txt"
write.table(summary_all_rounds, file = out_fn, sep = "\t", quote = FALSE, row.names = FALSE)

save.image(file = "conditional_test_outsideDRB1_"%&% gsub(",", "_", HLA_LABEL) %&%".RData")

