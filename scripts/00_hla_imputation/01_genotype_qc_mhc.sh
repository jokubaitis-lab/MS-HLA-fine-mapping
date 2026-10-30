#!/usr/bin/env bash
# Manuscript-oriented outline adapted from the study genotype QC command history.
# Requires pre-QC ancestry/relatedness decisions and reference harmonisation to be completed separately.
# Usage: bash 01_genotype_qc_mhc.sh INPUT_PLINK_PREFIX OUTPUT_DIR
set -euo pipefail
if [[ $# -ne 2 ]]; then echo "Usage: $0 INPUT_PLINK_PREFIX OUTPUT_DIR" >&2; exit 2; fi
input="$1"; out="$2"; mkdir -p "$out"
for tool in plink; do command -v "$tool" >/dev/null || { echo "Missing $tool" >&2; exit 1; }; done
# Generate initial per-sample and per-variant missingness summaries.
plink --bfile "$input" --missing --out "$out/initial_missingness"
# Initial permissive variant QC, matching the historical QC stage.
plink --bfile "$input" --geno 0.10 --hwe 1e-10 --make-bed --out "$out/initial_variant_qc"
# Final per-variant missingness and minor allele count filters.
# Apply to an input already cleaned for duplicates, strand, sample missingness,
# heterozygosity, relatedness, ancestry and reference allele frequency discrepancies.
# This script intentionally does NOT claim to perform those additional QC steps.
plink --bfile "$out/initial_variant_qc" --geno 0.02 --mac 5 --make-bed --out "$out/variant_qc"
# GRCh37 chr6:28-34 Mb was the HLA imputation preparation interval.
plink --bfile "$out/variant_qc" --chr 6 --from-mb 28 --to-mb 34 --make-bed --out "$out/hla_chr6_28_34Mb"
plink --bfile "$out/hla_chr6_28_34Mb" --freq --missing --out "$out/hla_chr6_28_34Mb_qc"
echo "Review and complete cohort-specific QC and reference harmonisation before imputation."
