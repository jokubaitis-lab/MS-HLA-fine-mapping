#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(stringr)
})

`%&%` <- function(a,b) paste0(a,b)

# ====================== Usage ======================
# Rscript flexible_haplo_multivar.R \
#   DICT_RDS \
#   RAW_PREFIX  ALLELE_NAME_TXT \
#   PHENO.txt PHENO_COL \
#   COVAR.txt  "PC1 + PC2 + PC3 + PC4 + PC5" \
#   TARGET_STR \
#   "DQB1:37+30;A:95;DPB1:65;B:116" \
#   OUT_PREFIX  [MIN_CARRIERS]
#
# TARGET_STR examples: "DRB1:71+13+86", "DQB1:37+30", "A:95"
# This script:
#   - Builds haplotypes for TARGET_STR (1+ positions).
#   - Reference haplotype = most frequent in controls.
#   - For each adjustor GENE:POS (from adjust_str), selects reference residue(s)
#     as most common in controls (ties allowed) and omits them from the model.
# ===================================================

args <- commandArgs(TRUE)
if (length(args) < 10) {
  stop("\nUsage:\n  Rscript flexible_haplo_multivar.R \\\n",
       "    DICT_RDS \\\n",
       "    RAW_PREFIX  ALLELE_NAME_TXT \\\n",
       "    PHENO.txt PHENO_COL \\\n",
       "    COVAR.txt  \"PC1 + PC2 + PC3 + PC4 + PC5\" \\\n",
       "    TARGET_STR \\\n",
       "    \"DQB1:37+30;A:95;DPB1:65;B:116\" \\\n",
       "    OUT_PREFIX  [MIN_CARRIERS]\n")
}

dict_rds     <- args[1]
raw_prefix   <- args[2]
allele_txt   <- args[3]
pheno_fp     <- args[4]
pheno_col    <- args[5]
covar_fp     <- args[6]
covar_terms  <- args[7]
target_str   <- args[8]   # NEW
adjust_str   <- args[9]
out_pref     <- args[10]
MIN_CARRIERS <- if (length(args) >= 11) as.integer(args[11]) else 0L

# 1-letter to 3-letter AA map (for labels)
aa3_map <- c(
  A="Ala", R="Arg", N="Asn", D="Asp", C="Cys",
  Q="Gln", E="Glu", G="Gly", H="His", I="Ile",
  L="Leu", K="Lys", M="Met", F="Phe", P="Pro",
  S="Ser", T="Thr", W="Trp", Y="Tyr", V="Val",
  `*`="*", X="X"
)

aa3 <- function(x) { x <- as.character(x); y <- aa3_map[x]; y[is.na(y)] <- x[is.na(y)]; unname(y) }

pretty_allele <- function(x) {
  x <- gsub("^HLA_[A-Z0-9]+_", "", x)
  x <- gsub("_", ":", x)
  paste0("*", x)
}

# ------------ Parse TARGET_STR and adjust_str -------------
parse_target <- function(s){
  # "DRB1:71+13+86" -> list(gene="DRB1", pos=c("71","13","86"))
  s <- trimws(s)
  gp <- unlist(strsplit(s, ":", fixed=TRUE))
  if (length(gp) != 2) stop("TARGET_STR must be like GENE:pos+pos+...")
  gene <- gp[1]
  pos  <- unlist(strsplit(gp[2], "\\+"))
  pos  <- trimws(pos[nzchar(pos)])
  if (!length(pos)) stop("No positions parsed from TARGET_STR.")
  list(gene = gene, pos = pos)
}

parse_adjust <- function(s){
  if (!nzchar(s)) return(list())
  parts <- trimws(unlist(strsplit(s, ";")))
  parts <- parts[nzchar(parts)]
  out <- lapply(parts, function(p){
    gp <- unlist(strsplit(p, ":", fixed=TRUE))
    if (length(gp) != 2) return(NULL)
    gene <- gp[1]
    poss <- trimws(unlist(strsplit(gp[2], "\\+")))
    list(gene=gene, pos=poss)
  })
  Filter(Negate(is.null), out)
}

TARGET <- parse_target(target_str)
TARGET_GENE <- TARGET$gene
TARGET_POS  <- TARGET$pos

adjust_spec <- parse_adjust(adjust_str)

# ------------------ Load resources -----------------
info_all <- readRDS(dict_rds)
need_cols <- c("gene","pos","AA","hla")
if (!all(need_cols %in% colnames(info_all))) {
  stop("Dictionary RDS must contain columns: ", paste(need_cols, collapse=", "))
}
info_all <- info_all %>%
  mutate(gene = as.character(gene),
         pos  = as.character(pos),
         AA   = as.character(AA),
         hla  = as.character(hla),
         tag  = paste0(gene, ":", pos, ":", AA))

# Dosage matrix
dose <- fread(raw_prefix %&% ".raw", header = TRUE, data.table = FALSE, check.names = FALSE)
if (ncol(dose) >= 6) dose <- dose[, c(2, 7:ncol(dose)), drop = FALSE]  # keep IID + markers
colnames(dose)[1] <- "IID"
allelenames <- read.table(allele_txt, header = FALSE, stringsAsFactors = FALSE)$V1
colnames(dose) <- c("IID", allelenames)

# Phenotype, covariates (assumes PCs named PC1..PC5; keeps what exists)
pheno <- read.table(pheno_fp, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
covar <- read.table(covar_fp, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
cov_keep <- unique(c("IID", intersect(colnames(covar), c("PC1","PC2","PC3","PC4","PC5","D_MSGP1_AU","D_MSGP2_AU","D_MSPGX_AU"))))
covar <- covar[, cov_keep, drop = FALSE]

# Merge
dat <- merge(dose, pheno, by = "IID")
dat <- merge(dat,  covar, by = "IID")

# Impute missing dosages to 0
allele_cols <- setdiff(colnames(dat), c("IID", colnames(pheno), colnames(covar)))
for (cn in allele_cols) if (is.numeric(dat[[cn]])) dat[[cn]][is.na(dat[[cn]])] <- 0

# Binary phenotype to 0/1
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

# ---------- Utilities ----------
# Map residue tag -> classical alleles
tag_to_alleles <- split(info_all$hla, info_all$tag)

add_residue_columns <- function(dat, tags, map_tag_to_alleles) {
  added <- character(0)
  for (tg in tags) {
    if (tg %in% colnames(dat)) { added <- c(added, tg); next }
    alset <- map_tag_to_alleles[[tg]]
    if (is.null(alset)) next
    present <- intersect(alset, colnames(dat))
    if (!length(present)) next
    dat[[tg]] <- if (length(present) == 1L) dat[[present]] else rowSums(dat[, present, drop = FALSE])
    added <- c(added, tg)
  }
  list(dat = dat, added = added)
}

choose_ref_tags_by_controls <- function(dat, tags, idx_ctrl) {
  avail <- intersect(tags, colnames(dat))
  if (!length(avail)) return(character(0))
  means_ctrl <- vapply(avail, function(cn) mean(dat[idx_ctrl, cn], na.rm = TRUE), numeric(1))
  if (all(is.na(means_ctrl))) return(character(0))
  m <- max(means_ctrl, na.rm = TRUE)
  names(which(means_ctrl == m))
}

# ===== 1) Build TARGET haplotypes (or single-position groups) =====
info_target <- info_all %>% filter(gene == TARGET_GENE, pos %in% TARGET_POS) %>%
  select(hla, pos, AA) %>% distinct()

if (!nrow(info_target)) stop("No rows in dictionary for TARGET_STR (", TARGET_GENE, ":", paste(TARGET_POS, collapse="+"), ").")

# Wide by position, keep alleles with all target positions defined
allele_pos_wide <- info_target %>%
  tidyr::pivot_wider(names_from = pos, values_from = AA) %>%
  mutate(across(all_of(TARGET_POS), as.character))

# Drop alleles missing any target position
for (p in TARGET_POS) allele_pos_wide <- allele_pos_wide %>% filter(!is.na(.data[[p]]))
if (!nrow(allele_pos_wide)) stop("No classical alleles with complete residues for TARGET_STR.")

# Build human-readable haplotype label and mapping
# Label uses 3-letter AAs joined by '-' in the order of TARGET_POS
res_cols <- lapply(TARGET_POS, function(p) aa3(allele_pos_wide[[p]]))
names(res_cols) <- TARGET_POS
haplo_label <- do.call(function(...) paste(..., sep = "-"), res_cols)

haplo_df <- cbind(hla = allele_pos_wide$hla, as.data.frame(res_cols, stringsAsFactors=FALSE))
haplo_df$haplo <- haplo_label

haplo_map <- split(haplo_df$hla, haplo_df$haplo)  # haplo label -> classical allele vector

# Column name builder e.g., HAP_DRB1_71_13_86_Ser-Arg-Val
pos_block <- paste(TARGET_POS, collapse = "_")
hap_colname <- function(h) paste0("HAP_", TARGET_GENE, "_", pos_block, "_", gsub("[^A-Za-z0-9]+", "_", h))

# Create haplotype dosage columns
dat_hap <- dat
hap_cols <- character(0)
for (h in names(haplo_map)) {
  alset <- unique(haplo_map[[h]])
  present <- intersect(alset, colnames(dat_hap))
  if (!length(present)) next
  cn <- hap_colname(h)
  dat_hap[[cn]] <- if (length(present) == 1L) dat_hap[[present]] else rowSums(dat_hap[, present, drop = FALSE])
  hap_cols <- c(hap_cols, cn)
}
if (!length(hap_cols)) {
  stop("No haplotype columns could be created for TARGET_STR. Check columns vs dictionary mapping.")
}

# Reference haplotype = most frequent in controls
af_ctrl_hap <- sapply(hap_cols, function(cn) mean(dat_hap[idx_ctrl, cn], na.rm = TRUE)/2)
ref_hap_col <- names(which.max(af_ctrl_hap))

# Optional rarity filter
if (MIN_CARRIERS > 0) {
  carriers <- sapply(hap_cols, function(cn) sum(dat_hap[[cn]] > 0, na.rm = TRUE))
  keep_haps <- names(which(carriers >= MIN_CARRIERS))
  drop_haps <- setdiff(hap_cols, keep_haps)
  if (length(drop_haps)) {
    message("Dropping rare haplotypes (< ", MIN_CARRIERS, " carriers): ", paste(drop_haps, collapse = ", "))
    hap_cols <- keep_haps
    if (!(ref_hap_col %in% hap_cols) && length(hap_cols) > 0) {
      af_ctrl_hap <- sapply(hap_cols, function(cn) mean(dat_hap[idx_ctrl, cn], na.rm = TRUE)/2)
      ref_hap_col <- names(which.max(af_ctrl_hap))
    }
  }
}
if (!length(hap_cols)) stop("All haplotypes dropped by MIN_CARRIERS.")

# ===== 2) Build adjusting AA residue terms for intended positions only =====
all_adjust_terms <- character(0)
ref_map <- list()   # (gene:pos) -> reference residue tag(s)

for (blk in adjust_spec) {
  adjust_gene <- blk$gene
  for (p in blk$pos) {
    pos_char <- as.character(p)
    tags <- info_all %>%
      filter(.data$gene == .env$adjust_gene, .data$pos == .env$pos_char) %>%
      transmute(tag = paste0(.data$gene, ":", .data$pos, ":", .data$AA)) %>%
      distinct() %>% pull(tag)
    if (!length(tags)) next

    tmp <- add_residue_columns(dat_hap, tags, tag_to_alleles)
    dat_hap <- tmp$dat
    avail <- intersect(tags, colnames(dat_hap))
    if (!length(avail)) next

    ref_tags <- choose_ref_tags_by_controls(dat_hap, avail, idx_ctrl)
    if (!length(ref_tags)) {
      means_all <- vapply(avail, function(cn) mean(dat_hap[, cn], na.rm = TRUE), numeric(1))
      m <- max(means_all, na.rm = TRUE)
      ref_tags <- names(which(means_all == m))
    }

    nonref <- setdiff(avail, ref_tags)
    ref_map[[paste0(gene,":",pos_char)]] <- ref_tags
    if (length(nonref)) all_adjust_terms <- c(all_adjust_terms, nonref)
  }
}

# Strict whitelist to ONLY intended positions
adj_pairs <- do.call(
  rbind,
  lapply(adjust_spec, function(blk) data.frame(gene=blk$gene, pos=blk$pos, stringsAsFactors=FALSE))
)
adj_pairs <- unique(adj_pairs)

keep_mask <- vapply(strsplit(all_adjust_terms, ":", fixed = TRUE), function(x) {
  if (length(x) != 3) return(FALSE)
  any(adj_pairs$gene == x[1] & adj_pairs$pos == x[2])
}, logical(1))
all_adjust_terms_strict <- unique(all_adjust_terms[keep_mask])

# ===== 3) Fit multivariate logistic model =====
hap_use <- setdiff(hap_cols, ref_hap_col)
if (!length(hap_use)) stop("No non-reference haplotypes available for modeling.")

dat_safe <- dat_hap
colnames(dat_safe) <- make.names(colnames(dat_safe))

adjust_terms_safe <- make.names(all_adjust_terms_strict)
hap_terms_safe    <- make.names(hap_use)

rhs_terms <- c(adjust_terms_safe, hap_terms_safe)
rhs <- paste(c(
  if (length(rhs_terms)) paste(rhs_terms, collapse = " + ") else NULL,
  covar_terms
), collapse = " + ")

cat("RHS:\n", rhs, "\n")

f <- as.formula(paste0(pheno_col, " ~ ", rhs))
print(f)
fit <- glm(f, data = dat_safe, family = binomial(link = "logit"))
co <- summary(fit)$coefficients

# ===== 4) Assemble output table (dynamic columns for positions) =====
# Map haplotype column -> human label and residues by position
hap_label_by_col <- setNames(names(haplo_map), sapply(names(haplo_map), hap_colname))  # col -> label like "Ser-Arg-Val"

split_label <- function(lab) strsplit(lab, "-", fixed = TRUE)[[1]]

# Helper to build a data.frame row for a given haplotype column
# --- replace your row_for_hap() with this base-R version (no !!!) ---
row_for_hap <- function(hcol, is_ref = FALSE) {
  lab <- hap_label_by_col[[hcol]]
  if (is.null(lab)) {
    lab <- gsub(paste0("^HAP_", TARGET_GENE, "_", pos_block, "_"), "", hcol)
    lab <- gsub("_", "-", lab)
  }
  parts <- strsplit(lab, "-", fixed = TRUE)[[1]]

  # AFs
  AF_ctrl <- mean(dat_hap[idx_ctrl, hcol], na.rm = TRUE) / 2
  AF_case <- mean(dat_hap[idx_case, hcol], na.rm = TRUE) / 2

  # OR/CI
  if (is_ref) {
    OR <- 1.0; L95 <- NA_real_; U95 <- NA_real_; pval <- NA_real_
  } else {
    safe <- make.names(hcol)
    if (!(safe %in% rownames(co))) return(NULL)
    beta <- co[safe, "Estimate"]; se <- co[safe, "Std. Error"]; pval <- co[safe, "Pr(>|z|)"]
    OR <- exp(beta); L95 <- exp(beta - 1.96 * se); U95 <- exp(beta + 1.96 * se)
  }

  # classical allele members
  al <- haplo_map[[lab]]
  al <- if (length(al)) paste(sort(pretty_allele(al)), collapse = ", ") else ""

  # Build a one-row data.frame with dynamic TARGET_POS columns
  pos_df <- as.data.frame(as.list(parts), stringsAsFactors = FALSE)
  names(pos_df) <- TARGET_POS

  # Bind stats columns
  stats_df <- data.frame(
    OR = OR, L95 = L95, U95 = U95, p = pval,
    AF_controls = AF_ctrl, AF_cases = AF_case,
    Classical_HLA_alleles = al,
    is_reference = is_ref,
    stringsAsFactors = FALSE
  )

  # cbind while preserving position column names
  out <- cbind(pos_df, stats_df, stringsAsFactors = FALSE)
  # keep raw names like "71", "13" (avoid X71)
  names(out)[seq_along(TARGET_POS)] <- TARGET_POS
  out
}

# --- and replace the block that builds out_tbl with this ---
ref_row <- row_for_hap(ref_hap_col, TRUE)
nonref_rows <- lapply(hap_use, function(h) row_for_hap(h, FALSE))
nonref_rows <- Filter(Negate(is.null), nonref_rows)

out_tbl <- dplyr::bind_rows(c(list(ref_row), nonref_rows)) %>%
  dplyr::arrange(is_reference, dplyr::desc(OR)) %>%
  dplyr::mutate(is_reference = ifelse(is_reference, "Reference", ""))

# Reorder columns: TARGET_POS first (in order), then stats
stat_cols <- c("OR","L95","U95","p","AF_controls","AF_cases","Classical_HLA_alleles","is_reference")
out_tbl <- out_tbl[, c(TARGET_POS, stat_cols), drop = FALSE]
# ===== 5) Write outputs with informative filenames =====
target_tag <- paste0(TARGET_GENE, "_", paste(TARGET_POS, collapse = "_"))
out_hap <- out_pref %&% ("_" %&% target_tag %&% "_multivar_table.tsv")
fwrite(out_tbl, out_hap, sep = "\t", quote = FALSE, na = "NA")
cat("Wrote haplotype table: ", out_hap, "\n", sep = "")

# Build readable reference maps for adjustor positions
ref_lines <- c()
if (length(ref_map)) {
  for (k in names(ref_map)) {
    refs <- sort(ref_map[[k]])
    if (!length(refs)) next
    residues <- paste0(gsub("^.+:.+:", "", refs), collapse = ",")
    ref_lines <- c(ref_lines, paste0("  ", k, " -> ref residue(s) in controls: ", residues,
                                     " [tags: ", paste(refs, collapse=", "), "]"))
  }
}

# Log details
out_log <- out_pref %&% ("_" %&% target_tag %&% "_model.txt")
writeLines(c(
  paste0("Multivariate logistic model for ", TARGET_GENE, " haplotype at positions: ", paste(TARGET_POS, collapse=", ")),
  paste0("Phenotype: ", pheno_col, " (0/1)"),
  paste0("Covariates: ", covar_terms),
  paste0("Target haplotype reference column: ", ref_hap_col),
  "Adjusting positions (non-reference residues ONLY):",
  if (length(all_adjust_terms_strict)) paste0("  ", paste(all_adjust_terms_strict, collapse = ", ")) else "  none",
  "Per-position adjustor reference residue(s) selected by control abundance:",
  if (length(ref_lines)) ref_lines else "  none",
  paste0("Formula: ", deparse(formula(fit)))
), con = out_log)
cat("Wrote model log: ", out_log, "\n", sep = "")

