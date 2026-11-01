#!/usr/bin/env bash
# Study-adapted HLA dosage processing; usage:
# bash 02_post_imputation_qc.sh CHR6_DOSE_VCF_GZ OUTPUT_PREFIX [SAMPLE_KEEP_FILE]
set -euo pipefail
[[ $# -ge 2 && $# -le 3 ]] || { echo "Usage: $0 VCF_GZ OUT [KEEP_FILE]" >&2; exit 2; }
vcf=$1; out=$2; keep=${3:-}
for tool in bcftools tabix plink2; do command -v "$tool" >/dev/null || exit 1; done
# Historical runs explored R2>0.5 and R2>0.7; 0.7 is the primary setting.
# Historical FILTER field contained PASS;GENOTYPED; check bcftools handling
# on the exact exported VCF rather than rewriting FILTER values blindly.
bcftools view -i 'INFO/R2>0.7' -Oz -o "${out}.r2_0.7.vcf.gz" "$vcf"
tabix -p vcf "${out}.r2_0.7.vcf.gz"
plink2 --vcf "${out}.r2_0.7.vcf.gz" dosage=DS --make-pgen --out "${out}.dosage"
# Do not guess FID/IID by splitting underscores; supply a reviewed mapping if needed.
plink2 --pfile "${out}.dosage" --missing variant-only --out "${out}.variant_missingness"
plink2 --pfile "${out}.dosage" --geno 0.05 --make-pgen --out "${out}.geno05"
plink2 --pfile "${out}.geno05" --missing sample-only --out "${out}.sample_missingness"
plink2 --pfile "${out}.geno05" --mind 0.05 --make-pgen --out "${out}.sample_qc"
if [[ -n "$keep" ]]; then
  [[ -f "$keep" ]] || { echo "Keep file not found" >&2; exit 1; }
  plink2 --pfile "${out}.sample_qc" --keep "$keep" --make-pgen --out "${out}.selected"
  selected="${out}.selected"
else selected="${out}.sample_qc"; fi
plink2 --pfile "$selected" --mac 5 --make-pgen --out "${out}.final_mac5"
plink2 --pfile "${out}.final_mac5" --export vcf bgz --out "${out}.final_mac5"
echo 'Confirm whether HLA alleles/amino acids are present in dosage VCF; their downstream encoding is not reproduced here.'
