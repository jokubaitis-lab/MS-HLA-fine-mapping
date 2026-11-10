#!/usr/bin/env bash
# Identify near-identical samples shared between ANZgene and MSBase cohorts.
# Requires PLINK 1.9; all inputs must already use the same genome build.
# Public adaptation of historical QC_filtered.sh; requires reviewed, harmonised inputs.
set -euo pipefail

usage() {
  cat <<'HELP'
Usage:
  bash QC_filtered.sh ANZGENE_PREFIX MSBASE_MANIFEST OUTPUT_DIR

Inputs:
  ANZGENE_PREFIX  PLINK binary prefix (.bed/.bim/.fam), already harmonised to GRCh37
  MSBASE_MANIFEST Tab-delimited file: cohort_name<TAB>plink_prefix (no header)
  OUTPUT_DIR      New or existing output directory

Example manifest:
  MSGP1_AU /path/to/MSGP1_AU
  MSGP1_CZ /path/to/MSGP1_CZ

Each comparison merges genotypes, removes duplicate variant positions,
excludes chr6:28-34 Mb from relatedness SNP selection,
prunes SNPs (MAF >= 0.05, --indep-pairwise 100 2 0.1), and calculates PI_HAT.
Only CROSS-COHORT pairs with PI_HAT > 0.98 are treated as duplicates.
The ANZgene individual is selected for removal; MSBase is retained.
Other pairs with PI_HAT > 0.1 are reported but not automatically removed.

Important: perform any required GRCh38-to-GRCh37 liftover upstream;
this script does not perform coordinate conversion. Verify genome build,
variant IDs, alleles and merge logs.
If PLINK reports unresolved merge conflicts, the script stops rather than
silently dropping variants or changing alleles.
HELP
}
[[ $# -eq 3 ]] || { usage >&2; exit 2; }
ANZ=$1; MANIFEST=$2; OUT=$3
command -v plink >/dev/null || { echo 'ERROR: PLINK 1.9 not found' >&2; exit 1; }
[[ -s "$MANIFEST" ]] || { echo 'ERROR: manifest not found' >&2; exit 1; }
for ext in bed bim fam; do [[ -s "${ANZ}.${ext}" ]] || { echo "ERROR: missing ${ANZ}.${ext}" >&2; exit 1; }; done
mkdir -p "$OUT" "$OUT/work" "$OUT/comparisons"
OUT=$(cd "$OUT" && pwd)
ANZ=$(cd "$(dirname "$ANZ")" && pwd)/$(basename "$ANZ")

# PLINK --update-ids needs OLD_FID OLD_IID NEW_FID NEW_IID.
# Preserve IID, replace FID with distinct cohort tag; preserve mappings for audit.
awk 'BEGIN{OFS="\t"} {print $1,$2,"ANZgene",$2}' "${ANZ}.fam" > "$OUT/work/anzgene_update_ids.tsv"
if [[ $(awk '{print $2}' "${ANZ}.fam" | sort | uniq -d | wc -l) -gt 0 ]]; then
  echo 'ERROR: ANZgene IIDs are not unique; cohort tagging by IID would be ambiguous.' >&2
  exit 1
fi
plink --bfile "$ANZ" --update-ids "$OUT/work/anzgene_update_ids.tsv" \
  --make-bed --out "$OUT/work/anzgene_tagged"

printf 'cohort\tFID_ANZgene\tIID_ANZgene\tFID_MSBase\tIID_MSBase\tPI_HAT\n' > "$OUT/cross_cohort_duplicates.tsv"
printf 'cohort\tFID1\tIID1\tFID2\tIID2\tPI_HAT\n' > "$OUT/related_pairs_gt_0.1.tsv"

while IFS=$'\t' read -r cohort msbase || [[ -n "${cohort:-}" ]]; do
  [[ -z "${cohort:-}" || "${cohort:0:1}" == '#' ]] && continue
  [[ "$cohort" =~ ^[A-Za-z0-9_-]+$ ]] || { echo "ERROR: invalid cohort label: $cohort" >&2; exit 1; }
  [[ -n "${msbase:-}" ]] || { echo "ERROR: missing prefix for $cohort" >&2; exit 1; }
  for ext in bed bim fam; do [[ -s "${msbase}.${ext}" ]] || { echo "ERROR: missing ${msbase}.${ext}" >&2; exit 1; }; done
  work="$OUT/work/$cohort"
  mkdir -p "$work"
  msbase=$(cd "$(dirname "$msbase")" && pwd)/$(basename "$msbase")
  if [[ $(awk '{print $2}' "${msbase}.fam" | sort | uniq -d | wc -l) -gt 0 ]]; then
    echo "ERROR: non-unique MSBase IIDs in $cohort" >&2
    exit 1
  fi
  awk -v c="MSBase_${cohort}" 'BEGIN{OFS="\t"} {print $1,$2,c,$2}' "${msbase}.fam" > "$work/msbase_update_ids.tsv"
  plink --bfile "$msbase" --update-ids "$work/msbase_update_ids.tsv" --make-bed --out "$work/msbase_tagged"

  # SNP-only filter, remove duplicate variant positions within each dataset.
  for spec in "anz:$OUT/work/anzgene_tagged" "ms:$work/msbase_tagged"; do
    tag=${spec%%:*}; source=${spec#*:}
    plink --bfile "$source" --snps-only just-acgt --make-bed --out "$work/${tag}_snps"
    plink --bfile "$work/${tag}_snps" --list-duplicate-vars suppress-first ids-only --out "$work/${tag}_dups"
    plink --bfile "$work/${tag}_snps" --exclude "$work/${tag}_dups.dupvar" --make-bed --out "$work/${tag}_nodup"
  done

  # Require a clean merge; unresolved allele/strand mismatches need manual QC.
  plink --bfile "$work/anz_nodup" \
    --bmerge "$work/ms_nodup.bed" "$work/ms_nodup.bim" "$work/ms_nodup.fam" \
    --make-bed --out "$work/merged" || {
      echo "ERROR: $cohort merge failed. Inspect $work/merged-merge.missnp and logs." >&2
      exit 1
    }
  # Exclude extended MHC from the IBD marker set to reduce long-range LD bias.
  # Retain the complete merged dataset for other downstream purposes.
  awk '$1==6 && $4>=28000000 && $4<=34000000 {print $2}' "$work/merged.bim" > "$work/mhc_exclude.snps"
  plink --bfile "$work/merged" --exclude "$work/mhc_exclude.snps" \
    --maf 0.05 --indep-pairwise 100 2 0.1 --out "$work/prune"
  plink --bfile "$work/merged" --extract "$work/prune.prune.in" --make-bed --out "$work/pruned"
  plink --bfile "$work/pruned" --genome --out "$work/ibd"

  # Header-based parsing avoids assumptions about column numbers.
  # Report only cross-cohort pairs, with the ANZgene individual listed first.
  awk -v c="$cohort" 'BEGIN{OFS="\t"}
    NR==1 {for(i=1;i<=NF;i++) h[$i]=i; next}
    {
      f1=$(h["FID1"]); f2=$(h["FID2"]); p=$(h["PI_HAT"])+0
      if (f1=="ANZgene" && f2=="MSBase_"c) {a1=f1;a2=$(h["IID1"]);b1=f2;b2=$(h["IID2"])}
      else if (f2=="ANZgene" && f1=="MSBase_"c) {a1=f2;a2=$(h["IID2"]);b1=f1;b2=$(h["IID1"])}
      else next
      if(p>0.98) print c,a1,a2,b1,b2,p
    }' "$work/ibd.genome" >> "$OUT/cross_cohort_duplicates.tsv"
  awk -v c="$cohort" 'BEGIN{OFS="\t"}
    NR==1 {for(i=1;i<=NF;i++) h[$i]=i; next}
    $(h["PI_HAT"])+0>0.1 {print c,$(h["FID1"]),$(h["IID1"]),$(h["FID2"]),$(h["IID2"]),$(h["PI_HAT"])}' \
    "$work/ibd.genome" >> "$OUT/related_pairs_gt_0.1.tsv"
done < "$MANIFEST"

# Convert tagged ANZgene IDs back to ORIGINAL FID/IID for PLINK --remove.
# Avoid silently removing all records when an IID is non-unique in ANZgene.
awk 'NR>1 {print $2,$3}' "$OUT/cross_cohort_duplicates.tsv" | sort -u > "$OUT/work/anzgene_tagged_remove.txt"
awk 'BEGIN{OFS="\t"} NR==FNR {remove[$1 FS $2]=1; next}
  (($3 FS $4) in remove) {print $1,$2}' \
  "$OUT/work/anzgene_tagged_remove.txt" "$OUT/work/anzgene_update_ids.tsv" \
  | sort -u > "$OUT/anzgene_remove_original_ids.txt"

# Produce a new ANZgene PLINK dataset, preserving original identifiers.
if [[ -s "$OUT/anzgene_remove_original_ids.txt" ]]; then
  plink --bfile "$ANZ" --remove "$OUT/anzgene_remove_original_ids.txt" \
    --make-bed --out "$OUT/ANZgene_deduplicated"
else
  plink --bfile "$ANZ" --make-bed --out "$OUT/ANZgene_deduplicated"
fi

echo "Completed. Review cross_cohort_duplicates.tsv and PLINK logs before downstream analysis."
echo "ANZgene removal list: $OUT/anzgene_remove_original_ids.txt"
