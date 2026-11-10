#!/usr/bin/env bash
# Stage 1: exploratory DRB1 allele selection; stage 2: primary MHC conditioning
# with the four verified DRB1 alleles stated in the manuscript.
# IMPORTANT: Selection and threshold conventions should be checked against final study logs.
set -euo pipefail
: "${PFILE:?Set PFILE}"
: "${PHENO_FILE:?Set PHENO_FILE}"
: "${COVAR_FILE:?Set COVAR_FILE}"
: "${INITIAL_GLM:?Set INITIAL_GLM to the single-marker .glm.logistic.hybrid output}"
PLINK2="${PLINK2:-plink2}"
PHENO_NAME="${PHENO_NAME:-MS_pheno}"
OUT_ROOT="${OUT_ROOT:-results/mhc_stepwise}"
COVAR_NAMES="${COVAR_NAMES:-PC1,PC2,PC3,PC4,PC5,D_MSGP1_AU,D_MSGP2_AU,D_MSPGX_AU}"
DRB1_THRESHOLD="${DRB1_THRESHOLD:-3.114100647732935e-06}"
MHC_THRESHOLD="${MHC_THRESHOLD:-3.114100647732935e-06}"
MAX_STEPS="${MAX_STEPS:-30}"
# For manuscript primary non-DRB1 analysis supply a reviewed four-allele list
# containing DRB1*15:01, *03:01, *13:03 and *01:03 in exact PLINK marker IDs.
# If omitted, the script stops before MHC-wide conditioning rather than silently
# treating an exploratory forward-selected DRB1 list as the published model.
DRB1_CONDITION_LIST="${DRB1_CONDITION_LIST:-}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$(dirname "$OUT_ROOT")"
# Select the smallest valid ADD P-value matching an ID prefix; emit ID and P.
select_top() {
  local file="$1" pattern="$2" threshold="$3" exclude="${4:-}"
  python3 - "$file" "$pattern" "$threshold" "$exclude" <<'PYSELECT'
import csv,math,re,sys
path,pattern,threshold,exclude=sys.argv[1:]
seen=set()
if exclude:
    with open(exclude) as f: seen={x.strip() for x in f if x.strip()}
best=None
with open(path,newline='') as f:
    r=csv.DictReader(f,delimiter='\t')
    if not r.fieldnames or not {'ID','P','TEST'}.issubset(r.fieldnames):
        sys.exit('Expected PLINK columns ID, P and TEST')
    for row in r:
        vid=row['ID']
        if row['TEST']!='ADD' or vid in seen or not re.search(pattern,vid):continue
        try:p=float(row['P'])
        except (TypeError,ValueError):continue
        if math.isfinite(p) and p<float(threshold) and (best is None or p<best[0]):best=(p,vid)
if best:print(best[1],best[0],sep='\t')
PYSELECT
}
run_conditioned() {
  local prefix="$1" conditions="$2"
  "$PLINK2" --pfile "$PFILE" --pheno "$PHENO_FILE" --pheno-name "$PHENO_NAME" \
    --covar "$COVAR_FILE" --covar-name "$COVAR_NAMES" \
    --glm omit-ref hide-covar cols=chrom,pos,ref,alt,test,nobs,beta,se,ci,tz,p,a1freqcc,a1freq \
    --ci 0.95 --covar-variance-standardize PC1 PC2 PC3 PC4 PC5 \
    --condition-list "$conditions" --out "$prefix"
  local result="${prefix}.${PHENO_NAME}.glm.logistic.hybrid"
  test -s "$result" || { echo "Missing PLINK result: $result" >&2; exit 1; }
  python3 "$SCRIPT_DIR/../shared/sort_plink_logistic.py" "$result" \
    "${prefix}.${PHENO_NAME}.glm.logistic.sorted.tsv"
}
# The original first stage selected HLA_DRB1 classical allele IDs.
DRB1_LIST="${OUT_ROOT}.condition_DRB1.list"
MHC_LIST="${OUT_ROOT}.condition_vars.list"
: > "$DRB1_LIST"
previous="$INITIAL_GLM"
for ((step=1;step<=MAX_STEPS;step++));do
  top="$(select_top "$previous" '^HLA[-_]?DRB1' "$DRB1_THRESHOLD" "$DRB1_LIST")"
  [[ -n "$top" ]] || break
  id="${top%%$'\t'*}"
  printf '%s\n' "$id" >> "$DRB1_LIST"
  prefix="${OUT_ROOT}.DRB1.step${step}"
  echo "DRB1 step ${step}: ${top}"
  run_conditioned "$prefix" "$DRB1_LIST"
  previous="${prefix}.${PHENO_NAME}.glm.logistic.hybrid"
done
# Primary MHC-wide conditioning uses the four manuscript-specified DRB1 alleles.
# Their exact variant IDs depend on the imputed PLINK dataset.
if [[ -z "$DRB1_CONDITION_LIST" || ! -s "$DRB1_CONDITION_LIST" ]]; then
  echo 'DRB1 forward-selection results were saved, but manuscript MHC conditioning requires DRB1_CONDITION_LIST with the four verified DRB1 allele IDs.' >&2
  exit 2
fi
mapfile -t reviewed_ids < <(grep -v '^[[:space:]]*#' "$DRB1_CONDITION_LIST" | awk 'NF {print $1}')
if [[ ${#reviewed_ids[@]} -ne 4 ]]; then
  echo 'DRB1_CONDITION_LIST must contain exactly four distinct manuscript DRB1 alleles.' >&2
  exit 2
fi
printf '%s\n' "${reviewed_ids[@]}" | sort -u > "$MHC_LIST"
if [[ $(wc -l < "$MHC_LIST") -ne 4 ]]; then
  echo 'DRB1_CONDITION_LIST has duplicate marker IDs.' >&2
  exit 2
fi
# Do not reuse the exploratory last-stage result: recompute the exact four-allele model.
run_conditioned "${OUT_ROOT}.DRB1_four_allele_baseline" "$MHC_LIST"
previous="${OUT_ROOT}.DRB1_four_allele_baseline.${PHENO_NAME}.glm.logistic.hybrid"
for ((step=1;step<=MAX_STEPS;step++));do
  top="$(select_top "$previous" '^(HLA_|AA_|rs)' "$MHC_THRESHOLD" "$MHC_LIST")"
  [[ -n "$top" ]] || break
  id="${top%%$'\t'*}"
  printf '%s\n' "$id" >> "$MHC_LIST"
  prefix="${OUT_ROOT}.step${step}"
  echo "MHC step ${step}: ${top}"
  run_conditioned "$prefix" "$MHC_LIST"
  previous="${prefix}.${PHENO_NAME}.glm.logistic.hybrid"
done
printf 'Final conditioning IDs:\n'; cat "$MHC_LIST"
