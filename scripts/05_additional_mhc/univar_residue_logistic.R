#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(stringr)
})

`%&%` <- function(a,b) paste0(a,b)
`%||%` <- function(x, y) if (is.null(x)) y else x   # <-- added

# ====================== Usage ======================
# Rscript univar_residue_logistic.R \
#   DICT_RDS \
#   RAW_PREFIX  ALLELE_NAME_TXT \
#   PHENO.txt PHENO_COL \
#   COVAR.txt  "PC1 + PC2 + PC3 + PC4 + PC5" \
#   TARGET_STR \
#   OUT_PREFIX  [COND_ALLELES]  [MIN_CARRIERS]
# ===================================================

args <- commandArgs(TRUE)
if (length(args) < 8) {
  stop("\nUsage:\n  Rscript univar_residue_logistic.R \\\n",
       "    DICT_RDS \\\n",
       "    RAW_PREFIX  ALLELE_NAME_TXT \\\n",
       "    PHENO.txt PHENO_COL \\\n",
       "    COVAR.txt  \"PC1 + PC2 + PC3 + PC4 + PC5 + batch\" \\\n",
       "    TARGET_STR \\\n",
       "    OUT_PREFIX  [COND_ALLELES]  [MIN_CARRIERS]\n")
}

dict_rds    <- args[1]
raw_prefix  <- args[2]
allele_txt  <- args[3]
pheno_fp    <- args[4]
pheno_col   <- args[5]
covar_fp    <- args[6]
covar_terms <- args[7]
target_str  <- args[8]
out_pref    <- args[9]
cond_alleles_str <- if (length(args) >= 10) args[10] else ""
MIN_CARRIERS <- if (length(args) >= 11) as.integer(args[11]) else 0L

aa3_map <- c(
  A="Ala", R="Arg", N="Asn", D="Asp", C="Cys",
  Q="Gln", E="Glu", G="Gly", H="His", I="Ile",
  L="Leu", K="Lys", M="Met", F="Phe", P="Pro",
  S="Ser", T="Thr", W="Trp", Y="Tyr", V="Val",
  `*`="*", X="X"
)
aa3 <- function(x){ x <- as.character(x); y <- aa3_map[x]; y[is.na(y)] <- x[is.na(y)]; unname(y) }

parse_target <- function(s){
  s <- trimws(s)
  gp <- unlist(strsplit(s, ":", fixed = TRUE))
  if (length(gp) != 2) stop("TARGET_STR must be like GENE:pos+pos+...")
  gene <- gp[1]
  pos  <- trimws(unlist(strsplit(gp[2], "\\+")))
  pos  <- pos[nzchar(pos)]
  if (!length(pos)) stop("No positions parsed from TARGET_STR.")
  list(gene = gene, pos = pos)
}

add_residue_columns <- function(dat, tags, tag_to_alleles) {
  added <- character(0)
  for (tg in tags) {
    if (tg %in% colnames(dat)) { added <- c(added, tg); next }
    alset <- tag_to_alleles[[tg]]
    if (is.null(alset)) next
    present <- intersect(alset, colnames(dat))
    if (!length(present)) next
    dat[[tg]] <- if (length(present) == 1L) dat[[present]] else rowSums(dat[, present, drop = FALSE])
    added <- c(added, tg)
  }
  list(dat = dat, added = added)
}

TARGET <- parse_target(target_str)
TARGET_GENE <- TARGET$gene
TARGET_POS  <- TARGET$pos

info_all <- readRDS(dict_rds)
need_cols <- c("gene","pos","AA","hla")
if (!all(need_cols %in% colnames(info_all))) {
  stop("Dictionary RDS must contain columns: gene,pos,AA,hla")
}
info_all <- info_all %>%
  mutate(gene = as.character(gene),
         pos  = as.character(pos),
         AA   = as.character(AA),
         hla  = as.character(hla),
         tag  = paste0(gene, ":", pos, ":", AA))

dose <- fread(raw_prefix %&% ".raw", header = TRUE, data.table = FALSE, check.names = FALSE)
if (ncol(dose) >= 6) dose <- dose[, c(2, 7:ncol(dose)), drop = FALSE]
colnames(dose)[1] <- "IID"
allelenames <- read.table(allele_txt, header = FALSE, stringsAsFactors = FALSE)$V1
colnames(dose) <- c("IID", allelenames)

pheno <- read.table(pheno_fp, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
covar <- read.table(covar_fp, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
cov_keep <- unique(c("IID", intersect(colnames(covar), c("PC1","PC2","PC3","PC4","PC5","D_MSGP1_AU","D_MSGP2_AU","D_MSPGX_AU"))))
covar <- covar[, cov_keep, drop = FALSE]

dat <- merge(dose, pheno, by = "IID")
dat <- merge(dat,  covar, by = "IID")

allele_cols <- setdiff(colnames(dat), c("IID", colnames(pheno), colnames(covar)))
for (cn in allele_cols) if (is.numeric(dat[[cn]])) dat[[cn]][is.na(dat[[cn]])] <- 0

if (is.numeric(dat[[pheno_col]]) &&
    max(dat[[pheno_col]], na.rm = TRUE) == 2 &&
    min(dat[[pheno_col]], na.rm = TRUE) == 1) {
  dat[[pheno_col]] <- dat[[pheno_col]] - 1
}
if (!all(dat[[pheno_col]] %in% c(0,1))) {
  stop("Phenotype column '", pheno_col, "' must be 0/1 or 1/2.")
}
idx_ctrl <- which(dat[[pheno_col]] == 0)
idx_case <- which(dat[[pheno_col]] == 1)

tag_to_alleles <- split(info_all$hla, info_all$tag)

target_residues <- info_all %>%
  filter(gene == TARGET_GENE, pos %in% TARGET_POS) %>%
  transmute(tag = paste0(gene, ":", pos, ":", AA),
            gene = gene, pos = pos, AA = AA) %>%
  distinct()
if (!nrow(target_residues)) {
  stop("No residues found for ", TARGET_GENE, " at positions ", paste(TARGET_POS, collapse=","))
}

tmp <- add_residue_columns(dat, target_residues$tag, tag_to_alleles)
dat <- tmp$dat
avail_tags <- intersect(target_residues$tag, colnames(dat))
if (!length(avail_tags)) stop("No residue dosage columns could be created (check dictionary vs .raw).")

if (MIN_CARRIERS > 0) {
  carriers <- sapply(avail_tags, function(cn) sum(dat[[cn]] > 0, na.rm = TRUE))
  keep <- names(which(carriers >= MIN_CARRIERS))
  if (!length(keep)) stop("All residues dropped by MIN_CARRIERS.")
  avail_tags <- keep
}

cond_alleles <- character(0)
if (nzchar(trimws(cond_alleles_str))) {
  cand <- trimws(unlist(strsplit(cond_alleles_str, ",")))
  cand <- cand[nzchar(cand)]
  cond_alleles <- intersect(cand, colnames(dat))
  missing <- setdiff(cand, cond_alleles)
  if (length(missing)) {
    message("Warning: conditioning allele columns not found and will be skipped: ",
            paste(missing, collapse = ", "))
  }
}

AF_controls <- sapply(avail_tags, function(cn) mean(dat[idx_ctrl, cn], na.rm = TRUE) / 2)
AF_cases    <- sapply(avail_tags, function(cn) mean(dat[idx_case, cn], na.rm = TRUE) / 2)

dat_safe <- dat
colnames(dat_safe) <- make.names(colnames(dat_safe))
pheno_safe <- make.names(pheno_col)
cond_safe <- make.names(cond_alleles)
rhs_cov <- trimws(covar_terms); rhs_cov <- if (nzchar(rhs_cov)) rhs_cov else NULL

fit_one_residue <- function(tag_raw){
  tag_safe <- make.names(tag_raw)
  rhs_parts <- c(tag_safe, cond_safe, rhs_cov)
  rhs <- paste(rhs_parts, collapse = " + ")
  f <- as.formula(paste0(pheno_safe, " ~ ", rhs))
  fit <- tryCatch(glm(f, data = dat_safe, family = binomial(link = "logit")),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  co <- suppressWarnings(summary(fit)$coefficients)
  if (!(tag_safe %in% rownames(co))) return(NULL)

  beta <- co[tag_safe, "Estimate"]
  se   <- co[tag_safe, "Std. Error"]
  pval <- co[tag_safe, "Pr(>|z|)"]
  OR   <- as.numeric(exp(beta))
  L95  <- as.numeric(exp(beta - 1.96*se))
  U95  <- as.numeric(exp(beta + 1.96*se))

  alset <- sort(unique(tag_to_alleles[[tag_raw]] %||% character(0)))
  pretty_allele <- function(x){
    x <- gsub("^HLA_", "", x); x <- gsub("_", ":", x); paste0("*", x)
  }
  al_str <- if (length(alset)) paste(pretty_allele(alset), collapse = ", ") else ""

  list(
    tag = tag_raw,
    OR = OR, L95 = L95, U95 = U95, p = pval,
    AF_controls = AF_controls[[tag_raw]],
    AF_cases    = AF_cases[[tag_raw]],
    Classical_HLA_alleles = al_str
  )
}

res_list <- lapply(avail_tags, fit_one_residue)
res_list <- Filter(Negate(is.null), res_list)
if (!length(res_list)) stop("All residue models failed (check predictors/covariates).")

out_df <- as.data.frame(do.call(rbind, lapply(res_list, as.data.frame)),
                        stringsAsFactors = FALSE)
tmp <- do.call(rbind, strsplit(out_df$tag, ":", fixed = TRUE))
colnames(tmp) <- c("Gene","Pos","AA")
out_df <- cbind(tmp, out_df[, setdiff(colnames(out_df), "tag"), drop = FALSE])
row.names(out_df) <- NULL

out_df <- out_df %>%
  dplyr::select(Gene, Pos, AA,
                OR, L95, U95, p,
                AF_controls, AF_cases, Classical_HLA_alleles)

tag_suffix <- paste0(TARGET_GENE, "_", paste(TARGET_POS, collapse = "_"))
out_tsv <- out_pref %&% ("_" %&% tag_suffix %&% "_univar_residues.tsv")
fwrite(out_df, out_tsv, sep = "\t", quote = FALSE, na = "NA")
cat("Wrote univariate residue results: ", out_tsv, "\n", sep = "")

out_log <- out_pref %&% ("_" %&% tag_suffix %&% "_univar_model.txt")
writeLines(c(
  paste0("Univariate residue logistic models for ", TARGET_GENE, " positions: ",
         paste(TARGET_POS, collapse = ", ")),
  paste0("Phenotype: ", pheno_col, " (0/1)"),
  paste0("Covariates (verbatim): ", if (is.null(rhs_cov)) "none" else rhs_cov),
  paste0("Conditioning classical allele columns: ",
         if (length(cond_alleles)) paste(cond_alleles, collapse=", ") else "none"),
  paste0("Residues modeled (N): ", length(avail_tags)),
  paste0("Carrier threshold MIN_CARRIERS: ", MIN_CARRIERS)
), con = out_log)
cat("Wrote model log: ", out_log, "\n", sep = "")

