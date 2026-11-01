#!/usr/bin/env bash
# Convert REVIEWED and reference-harmonised MHC PLINK dataset to upload VCF.
# Usage: bash 03_prepare_hla_imputation_input.sh INPUT_PLINK_PREFIX OUTPUT_PREFIX
set -euo pipefail
[[ $# -eq 2 ]] || { echo "Usage: $0 INPUT_PLINK_PREFIX OUTPUT_PREFIX" >&2; exit 2; }
for tool in plink bcftools tabix; do command -v "$tool" >/dev/null || exit 1; done
input=$1; out=$2
# IMPORTANT: original workflow used rename_bim.py and HRC-1000G-check-bim.pl
# with generated Run-plink.sh; these are reference-dependent, not reproducible
# without the exact panel and allele-check output. Run and review them upstream.
: "${REFERENCE_ALIGNMENT_CONFIRMED:?Set REFERENCE_ALIGNMENT_CONFIRMED=yes after checking build, IDs, strand and REF/ALT}"
[[ "$REFERENCE_ALIGNMENT_CONFIRMED" == yes ]] || exit 1
plink --bfile "$input" --recode vcf --out "$out"
bcftools sort "$out.vcf" -Oz -o "$out.vcf.gz"
tabix -p vcf "$out.vcf.gz"
echo 'Validate REF alleles against GRCh37 FASTA and inspect VCF before server submission.'
