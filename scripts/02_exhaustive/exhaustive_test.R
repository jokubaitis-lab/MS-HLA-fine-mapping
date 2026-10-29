#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

`%&%` <- function(a,b) paste0(a,b)

args <- commandArgs(TRUE)
if (length(args) < 10) {
  stop("\nUsage:\n  Rscript ra_exhaustive2.R \\\n",
       "    DICT_RDS \\\n",
       "    RAW_PREFIX  ALLELE_NAME_TXT \\\n",
       "    PHENO.txt PHENO_COL \\\n",
       "    COVAR.txt  \"PC1 + PC2 + PC3 + PC4 + PC5\" \\\n",
       "    GENE  Kmax  OUT_PREFIX  [ignore_positions]\n")
}

dict_rds   <- args[1]
raw_prefix <- args[2]
allele_txt <- args[3]
pheno_fp   <- args[4]
pheno_col  <- args[5]
covar_fp   <- args[6]
covars_str <- args[7]
gene       <- args[8]
Kmax       <- as.integer(args[9])
out_pref   <- args[10]
ignore_arg <- if (length(args) >= 11) args[11] else ""

# ---------------- Config ----------------
MIN_CARRIERS <- 10
# ----------------------------------------

# ---- Load dictionary & filter ----
info_all <- readRDS(dict_rds)
if (!all(c("gene","pos","AA","hla") %in% colnames(info_all))) {
  stop("Dictionary must contain columns: gene, pos, AA, hla")
}

info <- info_all %>%
  dplyr::filter(gene == !!gene) %>%
  dplyr::mutate(
    pos = trimws(as.character(pos)),
    tag = paste0(pos, ":", AA)
  )
if (!nrow(info)) stop("No dictionary rows for gene ", gene)

# ---- Robust ignore positions ----
if (nzchar(ignore_arg)) {
  toks <- unlist(strsplit(ignore_arg, "[,;\\s]+"))
  toks <- toks[nzchar(toks)]
  toks <- trimws(toks)

  expand_range <- function(x) {
    if (grepl("^-?\\d+\\s*-\\s*-?\\d+$", x)) {
      b <- as.integer(sub("\\s*-.*$", "", x))
      e <- as.integer(sub("^.*-\\s*", "", x))
      if (!is.na(b) && !is.na(e)) return(as.character(seq(min(b,e), max(b,e))))
    }
    x
  }
  toks <- unlist(lapply(toks, expand_range), use.names = FALSE)

  drop_pos <- unique(trimws(gsub("^.*?:", "", toks)))

  pos_char <- info$pos
  pos_num  <- suppressWarnings(as.integer(pos_char))
  drop_num <- suppressWarnings(as.integer(drop_pos))

  drop_hits <- pos_char %in% drop_pos
  if (any(!is.na(pos_num)) && any(!is.na(drop_num))) {
    drop_hits <- drop_hits | (pos_num %in% drop_num)
  }
  if (any(drop_hits)) {
    info <- info[!drop_hits, , drop=FALSE]
  }
  message("Ignoring ", gene, " positions: ", paste(unique(drop_pos), collapse=", "))
}

ALL_POS <- sort(unique(info$pos))
if (!length(ALL_POS)) stop("No positions left to test after ignoring.")
message("Remaining positions (n=", length(ALL_POS), "): ",
        paste(utils::head(ALL_POS, 20), collapse=", "),
        if (length(ALL_POS) > 20) ", ...")

# ---- Dosages / pheno / covars ----
dose <- fread(paste0(raw_prefix, ".raw"), check.names = FALSE)
if (!all(c("FID","IID") %in% names(dose))) {
  stop("Unexpected .raw format: missing FID/IID. Check input: ", raw_prefix, ".raw")
}
drop_cols <- intersect(c("FID","PAT","MAT","SEX","PHENOTYPE","PHENO"), names(dose))
dose <- dose[, setdiff(names(dose), drop_cols), with=FALSE]
setDF(dose)

allelenames <- scan(allele_txt, what="character", quiet=TRUE)
if (length(allelenames) != (ncol(dose)-1)) {
  stop("Allele-name count (", length(allelenames),
       ") does not match dosage columns (", ncol(dose)-1, ").")
}
colnames(dose) <- c("IID", allelenames)

pheno <- fread(pheno_fp) |> as.data.frame()
if (!("IID" %in% names(pheno))) stop("Phenotype file must contain IID column.")
if (!(pheno_col %in% names(pheno))) stop("Missing phenotype col: ", pheno_col)
if (is.numeric(pheno[[pheno_col]]) &&
    suppressWarnings(max(pheno[[pheno_col]],na.rm=TRUE)==2) &&
    suppressWarnings(min(pheno[[pheno_col]],na.rm=TRUE)==1)) {
  pheno[[pheno_col]] <- pheno[[pheno_col]] - 1
}
if (!all(na.omit(pheno[[pheno_col]]) %in% c(0,1))) {
  stop("Phenotype must be coded as 0/1 (or 1/2).")
}

covar <- fread(covar_fp) |> as.data.frame()
if (!("IID" %in% names(covar))) stop("Covariate file must contain IID column.")

# Align IIDs across tables
iid_common <- Reduce(intersect, list(dose$IID, pheno$IID, covar$IID))
if (!length(iid_common)) stop("No overlapping IID across dose/pheno/covar.")
dose  <- dose[dose$IID %in% iid_common, , drop=FALSE]
pheno <- pheno[pheno$IID %in% iid_common, , drop=FALSE]
covar <- covar[covar$IID %in% iid_common, , drop=FALSE]

# Covariates RHS (dedupe, keep only existing columns)
covars_str <- gsub("\\s+", "", covars_str)
cov_vec <- unique(strsplit(covars_str, "\\+")[[1]])
cov_vec <- cov_vec[nzchar(cov_vec)]
cov_vec <- intersect(cov_vec, colnames(covar))
covars_used <- if (length(cov_vec)) paste(cov_vec, collapse=" + ") else "1"

# Master merge of pheno + covars once
pheno_cov <- merge(pheno[, c("IID", pheno_col)],
                   covar[, c("IID", cov_vec), drop=FALSE],
                   by="IID", all=FALSE)
pheno_cov <- pheno_cov[complete.cases(pheno_cov), , drop=FALSE]

# Base: dose + pheno + covars (we’ll add group columns on copies)
base <- merge(dose, pheno_cov, by="IID", all=FALSE)

# ---- Grouping helpers ----
hla_by_tag <- split(info$hla, info$tag)

make_res_map <- function(info, pos_set) {
  positions <- unique(as.character(pos_set))
  if (!length(positions)) return(NULL)
  suppressWarnings(num <- as.numeric(positions))
  positions <- if (all(!is.na(num))) positions[order(num)] else sort(positions)
  tags_list <- lapply(positions, function(p) unique(info$tag[info$pos==p]))
  if (any(vapply(tags_list, length, integer(1))==0L)) return(NULL)
  comb <- do.call(expand.grid, c(tags_list, stringsAsFactors=FALSE))
  out <- vector("list", nrow(comb)); j <- 0L
  for (i in seq_len(nrow(comb))) {
    tags_vec <- as.character(unlist(comb[i, ], use.names=FALSE))
    inter <- Reduce(intersect, lapply(tags_vec, function(tg) hla_by_tag[[tg]]))
    if (!length(inter)) next
    j <- j + 1L
    out[[j]] <- data.frame(hap=paste(tags_vec, collapse="_"), hla=inter, stringsAsFactors=FALSE)
  }
  if (j==0L) return(NULL)
  do.call(rbind, out[seq_len(j)])
}

add_groups <- function(df, gm, alleles, min_car=MIN_CARRIERS) {
  adopted <- character(0)
  if (is.null(gm) || !nrow(gm)) return(list(dat=df, groups=adopted))
  gm <- gm[gm$hla %in% alleles, , drop=FALSE]
  if (!nrow(gm)) return(list(dat=df, groups=adopted))
  keep <- character(0)
  for (h in unique(gm$hap)) {
    pres <- intersect(gm$hla[gm$hap==h], colnames(df))
    if (!length(pres)) next
    carriers <- sum(rowSums(df[, pres, drop=FALSE]) > 0, na.rm=TRUE)
    if (carriers >= min_car) keep <- c(keep, h)
  }
  if (!length(keep)) return(list(dat=df, groups=adopted))
  for (h in keep) {
    pres <- intersect(gm$hla[gm$hap==h], colnames(df))
    df[[h]] <- if (length(pres)==1L) df[[pres]] else rowSums(df[, pres, drop=FALSE], na.rm=TRUE)
    adopted <- c(adopted, h)
  }
  list(dat=df, groups=adopted)
}

choose_ref <- function(df_with_groups, groups, pheno_cov, ph_col) {
  if (!length(groups)) return(NA_character_)
  tmp <- merge(df_with_groups[, c("IID", groups), drop=FALSE],
               pheno_cov[, c("IID", ph_col), drop=FALSE],
               by="IID", all.x=TRUE)
  ctrl <- subset(tmp, tmp[[ph_col]]==0); if (!nrow(ctrl)) ctrl <- tmp
  freqs <- vapply(groups, function(g) mean(ctrl[[g]] > 0, na.rm=TRUE), numeric(1))
  names(which.max(freqs))[1]
}

lrt_vs_null <- function(df_aug, groups, ph_col, covars_used, cov_vec) {
  need <- unique(c("IID", ph_col, cov_vec, groups))
  missing <- setdiff(need, colnames(df_aug))
  if (length(missing)) return(list(dev=NA, df=NA, p=NA))
  M <- df_aug[, need, drop=FALSE]
  M <- M[complete.cases(M), , drop=FALSE]
  if (!nrow(M)) return(list(dev=NA, df=NA, p=NA))
  colnames(M) <- make.names(colnames(M))
  ph_col_s <- make.names(ph_col)
  groups_s <- make.names(groups)
  has_cov  <- !(identical(covars_used, "1") || !length(cov_vec))
  rhs_groups <- if (length(groups_s)) paste(groups_s, collapse=" + ") else ""
  rhs <- if (has_cov && nzchar(rhs_groups)) {
    paste(rhs_groups, "+", covars_used)
  } else if (has_cov) {
    covars_used
  } else {
    rhs_groups
  }
  if (!nzchar(rhs)) return(list(dev=NA, df=NA, p=NA))
  f0 <- as.formula(paste0(ph_col_s, " ~ ", if (has_cov) covars_used else "1"))
  f1 <- as.formula(paste0(ph_col_s, " ~ ", rhs))
  suppressWarnings({
    m0 <- glm(f0, data=M, family=binomial())
    m1 <- glm(f1, data=M, family=binomial())
  })
  a <- tryCatch(anova(m0, m1, test="Chisq"), error=function(e) NULL)
  if (is.null(a) || nrow(a) < 2) return(list(dev=NA, df=NA, p=NA))
  list(dev=a$Deviance[2], df=a$Df[2], p=a$`Pr(>Chi)`[2])
}

# ---- NA-safe comparator for picking "best by Δdeviance" ----
better_than <- function(curr, best) {
  if (is.null(best)) return(TRUE)
  if (is.na(curr$dev) && is.na(best$dev)) return(FALSE)
  if (is.na(curr$dev)) return(FALSE)
  if (is.na(best$dev)) return(TRUE)
  curr$dev > best$dev
}

# ---- Exhaustive search over positions (WRITE PER-k) ----
best_rows <- list()
Kmax_eff <- max(1, min(Kmax, length(ALL_POS)))

for (k in 1:Kmax_eff) {
  message("Exhaustive: k = ", k)
  pos_combs <- utils::combn(ALL_POS, k, simplify=FALSE)

  rows_this_k <- vector("list", length(pos_combs))
  i <- 0L
  best <- NULL

  for (pcs in pos_combs) {
    gm  <- make_res_map(info, pcs)
    tmp <- add_groups(base, gm, colnames(dose)[-1], MIN_CARRIERS)

    if (!length(tmp$groups)) {
      i <- i + 1L
      rows_this_k[[i]] <- data.frame(
        gene=gene, k=k, combo=paste(pcs, collapse=","),
        dev=NA_real_, df=NA_real_, p=NA_real_,
        n_groups=0L, ref="", stringsAsFactors=FALSE
      )
      next
    }
    ref <- choose_ref(tmp$dat, tmp$groups, pheno_cov, pheno_col)
    use <- setdiff(tmp$groups, ref)
    if (!length(use)) {
      i <- i + 1L
      rows_this_k[[i]] <- data.frame(
        gene=gene, k=k, combo=paste(pcs, collapse=","),
        dev=NA_real_, df=NA_real_, p=NA_real_,
        n_groups=length(tmp$groups), ref=ifelse(is.na(ref),"",ref),
        stringsAsFactors=FALSE
      )
      next
    }

    st <- lrt_vs_null(tmp$dat, use, pheno_col, covars_used, cov_vec)
    row <- data.frame(
      gene  = gene,
      k     = k,
      combo = paste(pcs, collapse=","),
      dev   = st$dev,
      df    = st$df,
      p     = st$p,
      n_groups = length(use),
      ref   = ifelse(is.na(ref), "", ref),
      stringsAsFactors = FALSE
    )
    i <- i + 1L
    rows_this_k[[i]] <- row

    if (better_than(row, best)) best <- row
  }

  # bind and write per-k file (rank by Δdeviance desc; NA dev to bottom)
  tab_k <- if (i > 0) bind_rows(rows_this_k[seq_len(i)]) else
    data.frame(gene=character(), k=integer(), combo=character(),
               dev=numeric(), df=numeric(), p=numeric(),
               n_groups=integer(), ref=character(), stringsAsFactors=FALSE)

  if (nrow(tab_k)) {
    tab_k <- tab_k %>%
      mutate(
        dev_for_sort = ifelse(is.na(dev), -Inf, dev),
        p_for_sort   = ifelse(is.na(p),   Inf,  p)
      ) %>%
      arrange(desc(dev_for_sort), p_for_sort, desc(df)) %>%
      mutate(
        rank_by_delta = rank(-dev_for_sort, ties.method = "first"),
        rank_by_p     = rank(p_for_sort,    ties.method = "first")
      ) %>%
      select(-dev_for_sort, -p_for_sort)
  }

  out_k <- out_pref %&% "_k" %&% k %&% ".tsv"
  fwrite(tab_k, out_k, sep="\t", na="NA", quote=FALSE)
  message("Wrote per-k table (ranked by Δdeviance): ", out_k)

  if (!is.null(best)) {
    best_rows[[as.character(k)]] <- best
    message("  best k=", k, " combo=", best$combo,
            "  Δdev=", signif(best$dev, 4),
            "  p=", ifelse(is.na(best$p), "NA", signif(best$p, 3)))
  } else {
    message("  no valid combos for k=", k, " (likely no carriers or singular fits)")
  }
}

# Best-per-k summary
out_tab <- if (length(best_rows)) {
  data.table::rbindlist(best_rows, use.names=TRUE, fill=TRUE)
} else {
  data.frame(gene=character(), k=integer(), combo=character(),
             dev=numeric(), df=numeric(), p=numeric(),
             n_groups=integer(), ref=character(), stringsAsFactors = FALSE)
}
outfile <- out_pref %&% "_exhaustive_best.tsv"
fwrite(out_tab, outfile, sep="\t", na="NA", quote=FALSE)
message("Wrote best-per-k summary: ", outfile)

