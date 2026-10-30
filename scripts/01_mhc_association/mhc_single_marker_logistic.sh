#!/usr/bin/env bash
# Adapted from HLA_association_binary.cmd, Single-marker test.
set -euo pipefail
: "${PFILE:?Set PFILE to the QCed imputed MHC PLINK2 prefix}"
: "${PHENO_FILE:?Set PHENO_FILE}"
: "${COVAR_FILE:?Set COVAR_FILE}"
PLINK2="${PLINK2:-plink2}"
PHENO_NAME="${PHENO_NAME:-MS_pheno}"
OUT_PREFIX="${OUT_PREFIX:-results/mhc_single_marker}"
COVAR_NAMES="${COVAR_NAMES:-PC1,PC2,PC3,PC4,PC5,D_MSGP1_AU,D_MSGP2_AU,D_MSPGX_AU}"
mkdir -p "$(dirname "$OUT_PREFIX")"
"$PLINK2" --pfile "$PFILE" --pheno "$PHENO_FILE" --pheno-name "$PHENO_NAME" \
  --covar "$COVAR_FILE" --covar-name "$COVAR_NAMES" \
  --glm omit-ref hide-covar cols=chrom,pos,ref,alt,test,nobs,beta,se,ci,tz,p,a1freqcc,a1freq \
  --ci 0.95 --covar-variance-standardize PC1 PC2 PC3 PC4 PC5 --out "$OUT_PREFIX"
# PLINK2 output schema varies with version. Sort ADD rows by the named P column.
python3 "$(dirname "$0")/../shared/sort_plink_logistic.py" \
  "${OUT_PREFIX}.${PHENO_NAME}.glm.logistic.hybrid" \
  "${OUT_PREFIX}.${PHENO_NAME}.glm.logistic.sorted.tsv"
