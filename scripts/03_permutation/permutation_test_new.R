#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(stringr)
})

`%&%` <- function(a,b) paste0(a,b)

args <- commandArgs(TRUE)
if (length(args) < 12) {
  stop("\nUsage:\n  Rscript ra_permutation.R \\\n",
       "    DICT_RDS \\\n",
       "    RAW_PREFIX  ALLELE_NAME_TXT \\\n",
       "    PHENO.txt PHENO_COL \\\n",
       "    COVAR.txt  \"PC1 + PC2 + PC3 + PC4 + PC5\" \\\n",
       "    GENE  Kmax  Nperm  OUT_PREFIX \\\n",
       "    [ignore] [mode]\n\n",
       "  ignore: positions and/or residue tags separated by comma/space/semicolon; examples:\n",
       "          \"96 71\"   or  \"DRB1:96,DRB1:71\"   (drop whole positions)\n",
       "          \"71:K\"    or  \"DRB1:71:K\"        (drop specific residue tag)\n",
       "  mode:   sequential(default) | exhaustive\n")
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
Nperm      <- as.integer(args[10])
out_pref   <- args[11]
ignore_arg <- if (length(args) >= 12) args[12] else ""
SEL_MODE   <- if (length(args) >= 13) tolower(args[13]) else "sequential"
if (!SEL_MODE %in% c("exhaustive","sequential"))
  stop("mode must be 'exhaustive' or 'sequential'")

MIN_CARRIERS <- 10

# ---------------- Dictionary & positions ----------------
info_all <- readRDS(dict_rds)
info <- info_all %>% filter(.data$gene == !!gene) %>%
  mutate(pos = as.character(pos), tag = paste0(pos, ":", AA))
if (!nrow(info)) stop("No dict rows for ", gene)

# Ignore filter: positions and/or residue tags
if (nzchar(ignore_arg)) {
  toks <- unlist(strsplit(ignore_arg, "[,;\\s]+"))
  toks <- toks[nzchar(toks)]
  is_pos <- grepl("^([A-Za-z0-9_]+:)?\\d+$", toks)
  is_tag <- grepl("^([A-Za-z0-9_]+:)?\\d+:[A-Za-z]+$", toks)
  drop_pos <- unique(gsub("^.*?:", "", toks[is_pos]))                   # "96"
  drop_tag <- unique(gsub("^[A-Za-z0-9_]+:", "", toks[is_tag]))         # "71:K"
  if (length(drop_pos)) {
    info <- info[!(info$pos %in% drop_pos), , drop=FALSE]
    message("Ignoring ", gene, " positions: ", paste(drop_pos, collapse=", "))
  }
  if (length(drop_tag)) {
    info <- info[!(info$tag %in% drop_tag), , drop=FALSE]
    message("Ignoring ", gene, " residue tags: ", paste(drop_tag, collapse=", "))
  }
}
ALL_POS <- sort(unique(info$pos))
if (!length(ALL_POS)) stop("No positions left to test after ignore filter.")

# ---------------- Dosage / pheno / covar / base ----------------
dose <- fread(paste0(raw_prefix, ".raw"), check.names=FALSE)
# .raw typical columns: FID IID PAT MAT SEX <allele...>
dose <- dose[, -c(1,3:6), with=FALSE]  # keep IID and allele dosage columns
allelenames <- scan(allele_txt, what="character", quiet=TRUE)
setnames(dose, c("IID", allelenames))

pheno <- fread(pheno_fp)
if (!(pheno_col %in% names(pheno))) stop("Missing phenotype col: ", pheno_col)
# convert 1/2 -> 0/1 if needed
if (is.numeric(pheno[[pheno_col]]) && max(pheno[[pheno_col]],na.rm=TRUE)==2 &&
    min(pheno[[pheno_col]],na.rm=TRUE)==1) pheno[[pheno_col]] <- pheno[[pheno_col]] - 1

covar <- fread(covar_fp)
if (!("IID" %in% names(covar))) stop("Covariates need IID.")

iid_common <- Reduce(intersect, list(dose$IID, pheno$IID, covar$IID))
dose  <- as.data.frame(dose[dose$IID %in% iid_common, , drop=FALSE])
pheno <- as.data.frame(pheno[pheno$IID %in% iid_common, , drop=FALSE])
covar <- as.data.frame(covar[covar$IID %in% iid_common, , drop=FALSE])

covars_str <- gsub("\\s+","", covars_str)
cov_vec <- unique(strsplit(covars_str, "\\+")[[1]])
cov_vec <- cov_vec[nzchar(cov_vec)]
cov_vec <- intersect(cov_vec, colnames(covar))
covars_used <- if (length(cov_vec)) paste(cov_vec, collapse=" + ") else "1"

pheno_cov <- merge(pheno[, c("IID", pheno_col)], covar[, c("IID", cov_vec), drop=FALSE],
                   by="IID", all=FALSE)
pheno_cov <- pheno_cov[complete.cases(pheno_cov), , drop=FALSE]
base <- merge(as.data.frame(dose), pheno_cov, by="IID", all=FALSE)

# ---------------- Shared helpers ----------------
add_groups <- function(df, gm, alleles, min_car=MIN_CARRIERS) {
  adopted <- character(0)
  if (is.null(gm) || !nrow(gm)) return(list(dat=df, groups=adopted))
  gm <- gm[gm$hla %in% alleles, , drop=FALSE]
  if (!nrow(gm)) return(list(dat=df, groups=adopted))
  keep <- character(0)
  for (h in unique(gm$hap)) {
    pres <- intersect(gm$hla[gm$hap==h], colnames(df))
    if (!length(pres)) next
    carriers <- sum(rowSums(df[, pres, drop=FALSE])>0, na.rm=TRUE)
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
               pheno_cov[, c("IID", ph_col)], by="IID", all.x=TRUE)
  ctrl <- subset(tmp, tmp[[ph_col]]==0)
  if (!nrow(ctrl)) ctrl <- tmp
  freqs <- vapply(groups, function(g) mean(ctrl[[g]]>0, na.rm=TRUE), numeric(1))
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
  rhs_groups <- if (length(groups)) paste(make.names(groups), collapse=" + ") else NULL
  rhs <- paste(c(rhs_groups, covars_used), collapse=" + ")
  rhs <- gsub("^ \\+ | \\+ $", "", rhs)
  f0  <- as.formula(paste0(make.names(ph_col), " ~ ", covars_used))
  f1  <- as.formula(paste0(make.names(ph_col), " ~ ", rhs))
  m0 <- glm(f0, data=M, family=binomial())
  m1 <- glm(f1, data=M, family=binomial())
  a  <- anova(m0, m1, test="Chisq")
  list(dev=a$Deviance[2], df=a$Df[2], p=a$`Pr(>Chi)`[2])
}

make_res_map <- function(hla_by_tag, positions) {
  positions <- unique(as.character(positions))
  suppressWarnings(num <- as.numeric(positions))
  positions <- if (all(!is.na(num))) positions[order(num)] else sort(positions)
  tags_list <- lapply(positions, function(p) {
    grep(paste0("^", p, ":"), names(hla_by_tag), value=TRUE)
  })
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

# Compare-rows helper: best by largest deviance (ties: smaller df, then smaller p, then combo)
better_than <- function(a, b) {
  if (is.null(b)) return(TRUE)
  if (any(!(c("dev","df","p","combo") %in% names(a)))) return(FALSE)
  if (any(!(c("dev","df","p","combo") %in% names(b)))) return(TRUE)
  if (is.na(a$dev) && !is.na(b$dev)) return(FALSE)
  if (!is.na(a$dev) && is.na(b$dev)) return(TRUE)
  if (is.na(a$dev) && is.na(b$dev))  return(FALSE)
  if (a$dev != b$dev) return(a$dev > b$dev)
  if (!is.na(a$df) && !is.na(b$df) && a$df != b$df) return(a$df < b$df)
  if (!is.na(a$p)  && !is.na(b$p)  && a$p  != b$p)  return(a$p  < b$p)
  return(as.character(a$combo) < as.character(b$combo))
}

best_by_k_exhaustive <- function(hla_by_tag, ALL_POS, base, allelenames, pheno_cov, pheno_col, covars_used, cov_vec, Kmax) {
  out <- list()
  for (k in 1:min(Kmax, length(ALL_POS))) {
    best <- NULL
    pos_combs <- utils::combn(ALL_POS, k, simplify=FALSE)
    for (pcs in pos_combs) {
      gm <- make_res_map(hla_by_tag, pcs)
      tmp <- add_groups(base, gm, allelenames, MIN_CARRIERS)
      if (!length(tmp$groups)) next
      ref <- choose_ref(tmp$dat, tmp$groups, pheno_cov, pheno_col)
      use <- setdiff(tmp$groups, ref)
      if (!length(use)) next
      st <- lrt_vs_null(tmp$dat, use, pheno_col, covars_used, cov_vec)
      row <- data.frame(k=k, combo=paste(pcs, collapse=","), dev=st$dev, df=st$df, p=st$p)
      if (better_than(row, best)) best <- row
    }
    if (!is.null(best)) out[[k]] <- best
  }
  if (length(out)) do.call(rbind, out) else data.frame()
}

best_by_k_sequential <- function(hla_by_tag, ALL_POS, base, allelenames, pheno_cov, pheno_col, covars_used, cov_vec, Kmax) {
  chosen <- character(0)
  out <- list()
  for (k in 1:min(Kmax, length(ALL_POS))) {
    cand <- setdiff(ALL_POS, chosen)
    best <- NULL; best_pos <- NA_character_
    for (p in cand) {
      pcs <- c(chosen, p)
      gm  <- make_res_map(hla_by_tag, pcs)
      tmp <- add_groups(base, gm, allelenames, MIN_CARRIERS)
      if (!length(tmp$groups)) next
      ref <- choose_ref(tmp$dat, tmp$groups, pheno_cov, pheno_col)
      use <- setdiff(tmp$groups, ref)
      if (!length(use)) next
      st <- lrt_vs_null(tmp$dat, use, pheno_col, covars_used, cov_vec)
      row <- data.frame(k=k, combo=paste(pcs, collapse=","), dev=st$dev, df=st$df, p=st$p)
      if (better_than(row, best)) { best <- row; best_pos <- p }
    }
    if (is.null(best)) break
    chosen <- c(chosen, best_pos)
    out[[k]] <- best
  }
  if (length(out)) do.call(rbind, out) else data.frame()
}

all_combos_by_k <- function(hla_by_tag, ALL_POS, base, allelenames, pheno_cov, pheno_col, covars_used, cov_vec, Kmax) {
  out <- list()
  for (k in 1:min(Kmax, length(ALL_POS))) {
    rows <- list()
    pos_combs <- utils::combn(ALL_POS, k, simplify=FALSE)
    for (pcs in pos_combs) {
      gm  <- make_res_map(hla_by_tag, pcs)
      tmp <- add_groups(base, gm, allelenames, MIN_CARRIERS)
      if (!length(tmp$groups)) next
      ref <- choose_ref(tmp$dat, tmp$groups, pheno_cov, pheno_col)
      use <- setdiff(tmp$groups, ref)
      if (!length(use)) next
      st  <- lrt_vs_null(tmp$dat, use, pheno_col, covars_used, cov_vec)
      rows[[length(rows)+1L]] <- data.frame(
        k = k, combo = paste(pcs, collapse=","), dev = st$dev, df = st$df, p = st$p
      )
    }
    if (length(rows)) out[[length(out)+1L]] <- do.call(rbind, rows)
  }
  if (length(out)) do.call(rbind, out) else data.frame()
}

# ---------------- Observed: best and all ----------------
hla_by_tag_obs <- split(info$hla, info$tag)

obs_best <- if (SEL_MODE == "sequential") {
  best_by_k_sequential(hla_by_tag_obs, ALL_POS, base, allelenames, pheno_cov, pheno_col, covars_used, cov_vec, Kmax)
} else {
  best_by_k_exhaustive(hla_by_tag_obs, ALL_POS, base, allelenames, pheno_cov, pheno_col, covars_used, cov_vec, Kmax)
}
fwrite(obs_best, out_pref %&% "_perm_obs_best.tsv", sep="\t")
message("Wrote: ", out_pref %&% "_perm_obs_best.tsv")

obs_all <- all_combos_by_k(hla_by_tag_obs, ALL_POS, base, allelenames, pheno_cov, pheno_col, covars_used, cov_vec, Kmax)
if (nrow(obs_all)) {
  fwrite(obs_all, out_pref %&% "_perm_obs_all.tsv", sep="\t")
  for (kk in sort(unique(obs_all$k))) {
    fwrite(subset(obs_all, k==kk), out_pref %&% "_perm_obs_all_k" %&% kk %&% ".tsv", sep="\t")
  }
  message("Wrote: ", out_pref %&% "_perm_obs_all*.tsv")
}

# ---------------- Prepare allele x position (for permutations) ----------------
# NOTE: Built from 'info' AFTER ignore filter, so ignored positions/tags are excluded in permutations too.
allele_wide <- info %>%
  select(hla, pos, AA) %>% distinct() %>%
  tidyr::pivot_wider(names_from = pos, values_from = AA)

# ---------------- Permutations ----------------
perm_counts <- setNames(rep(0L, nrow(obs_best)), paste0("k", obs_best$k))
# store best dev/df/p per permutation for each k
perm_beststats <- list()  # list of data.frames per k: columns perm, k, dev_best, df_best, p_best

set.seed(1)
for (b in 1:Nperm) {
  perm_ix <- sample.int(nrow(allele_wide))
  perm_map <- allele_wide
  perm_map$hla <- allele_wide$hla
  # Shuffle each position independently, respecting ignore-filtered universe
  perm_map[, -1] <- allele_wide[perm_ix, -1, drop=FALSE]

  info_perm <- perm_map %>%
    tidyr::pivot_longer(-hla, names_to="pos", values_to="AA") %>%
    filter(!is.na(AA)) %>%
    mutate(tag = paste0(pos, ":", AA))

  hla_by_tag_perm <- split(info_perm$hla, info_perm$tag)

  pbest <- if (SEL_MODE == "sequential") {
    best_by_k_sequential(hla_by_tag_perm, ALL_POS, base, allelenames, pheno_cov, pheno_col, covars_used, cov_vec, Kmax)
  } else {
    best_by_k_exhaustive(hla_by_tag_perm, ALL_POS, base, allelenames, pheno_cov, pheno_col, covars_used, cov_vec, Kmax)
  }

  if (nrow(pbest)) {
    # record best dev/df/p per k (for density / summary plots)
    for (kk in unique(pbest$k)) {
      rowkk <- pbest[pbest$k == kk, , drop=FALSE]
      if (nrow(rowkk) && !is.na(rowkk$dev[1])) {
        key <- paste0("k", kk)
        add <- data.frame(
          perm = b,
          k = kk,
          dev_best = rowkk$dev[1],
          df_best  = rowkk$df[1],
          p_best   = rowkk$p[1]
        )
        if (is.null(perm_beststats[[key]])) {
          perm_beststats[[key]] <- add
        } else {
          perm_beststats[[key]] <- rbind(perm_beststats[[key]], add)
        }
      }
    }
    # count exceedances vs observed best (compare deviance only)
    for (i in seq_len(nrow(obs_best))) {
      k  <- obs_best$k[i]; dev_obs <- obs_best$dev[i]
      rowp <- pbest[pbest$k==k, , drop=FALSE]
      if (nrow(rowp) && !is.na(rowp$dev[1]) && !is.na(dev_obs) && rowp$dev[1] >= dev_obs) {
        perm_counts[[paste0("k", k)]] <- perm_counts[[paste0("k", k)]] + 1L
      }
    }
  }
  if (b %% 100 == 0) message("...perm ", b, "/", Nperm)
}

# Save permutation best-stat distributions (now includes p alongside dev and df)
if (length(perm_beststats)) {
  for (nm in names(perm_beststats)) {
    dfk <- perm_beststats[[nm]]
    if (nrow(dfk)) {
      # Keep original filename pattern but now with dev, df, and p columns
      fwrite(dfk, out_pref %&% "_perm_bestdev_" %&% nm %&% ".tsv", sep="\t")
    }
  }
  message("Wrote: ", out_pref %&% "_perm_bestdev_k*.tsv (columns: perm,k,dev_best,df_best,p_best)")
}

# Empirical p for observed-best per k (by deviance)
emp <- data.frame(
  k = obs_best$k,
  combo_obs = obs_best$combo,
  dev_obs = obs_best$dev,
  df_obs  = obs_best$df,
  p_obs   = obs_best$p,
  exceed  = as.integer(perm_counts[paste0("k", obs_best$k)]),
  nperm   = Nperm,
  p_emp   = (as.integer(perm_counts[paste0("k", obs_best$k)]) + 1) / (Nperm + 1)
)
fwrite(emp, out_pref %&% "_perm_empirical.tsv", sep="\t")
message("Wrote: ", out_pref %&% "_perm_empirical.tsv")

