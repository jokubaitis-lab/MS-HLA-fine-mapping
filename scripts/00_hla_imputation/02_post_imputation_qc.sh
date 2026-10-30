#!/usr/bin/env bash
# Study-oriented post-HLA-imputation processing of a Michigan Imputation Server dosage VCF.
# Usage: bash 02_post_imputation_qc.sh CHR6_DOSE_VCF_GZ OUTPUT_PREFIX
# Input is a GRCh37 HLA-imputed VCF containing INFO/R2 and genotype dosages.
set -euo pipefail
if [[ $# -ne 2 ]]; then echo "Usage: $0 CHR6_DOSE_VCF_GZ OUTPUT_PREFIX" >&2; exit 2; fi
vcf="$1"; out="$2"
for tool in bcftools tabix plink2; do command -v "$tool" >/dev/null || { echo "Missing $tool" >&2; exit 1; }; done
# Filter on imputation R2 > 0.7, the manuscript threshold.
# The historical script normalised PASS;GENOTYPED in FILTER before bcftools filtering.
bcftools view -i 'INFO/R2>0.7' -Oz -o "${out}.r2_0.7.vcf.gz" "$vcf"
tabix -p vcf "${out}.r2_0.7.vcf.gz"
plink2 --vcf "${out}.r2_0.7.vcf.gz" dosage=DS --make-pgen --out "${out}.dosage"
# MAC >= 5 as in the manuscript; additional sample/variant QC is cohort-specific.
plink2 --pfile "${out}.dosage" --mac 5 --make-pgen --out "${out}.mac5"
plink2 --pfile "${out}.mac5" --export vcf bgz --out "${out}.mac5"
echo "Downstream HLA allele/amino-acid dosage exports and sample-ID harmonisation are study-specific."
