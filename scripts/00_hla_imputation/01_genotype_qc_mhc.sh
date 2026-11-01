#!/usr/bin/env bash
# Reconstructed HLA-region genotype QC from historical study commands.
# Usage: bash 01_genotype_qc_mhc.sh INPUT_PLINK_PREFIX OUTPUT_DIR REFERENCE_FASTA
# Requires: PLINK 1.9, SNPflip, Python3, and project helper scripts.
# Important: sample relatedness/heterozygosity and frequency checks are REVIEW GATES.
set -euo pipefail
[[ $# -eq 3 ]] || { echo "Usage: $0 INPUT_PLINK_PREFIX OUTPUT_DIR REFERENCE_FASTA" >&2; exit 2; }
input=$1; out=$2; fasta=$3
mkdir -p "$out"
for tool in plink snpflip python3; do command -v "$tool" >/dev/null || { echo "Missing $tool" >&2; exit 1; }; done
: "${DUP_SCRIPT:?Set DUP_SCRIPT to reviewed get_duprem_var.py path}"
: "${HIGH_HET_SCRIPT:?Set HIGH_HET_SCRIPT to reviewed HighHeterozygosityOutlier.R path}"
: "${R_SCRIPT:=Rscript}"
plink --bfile "$input" --missing --out "$out/initial"
# Helper reads PLINK missingness and produces prefix.remdup.snp
python3 "$DUP_SCRIPT" "$out/initial"
[[ -f "$out/initial.remdup.snp" ]] || { echo 'Duplicate exclusion file missing' >&2; exit 1; }
plink --bfile "$input" --exclude "$out/initial.remdup.snp" --make-bed --out "$out/dedup"
snpflip --fasta-genome "$fasta" --bim-file "$out/dedup.bim" -o "$out/strand"
# Historical commands excluded ambiguous variants; flagged strand flips require
# separate correction/validation before proceeding, NOT automatic exclusion.
[[ -f "$out/strand.ambiguous" ]] || { echo 'SNPflip ambiguous list missing' >&2; exit 1; }
plink --bfile "$out/dedup" --exclude "$out/strand.ambiguous" --make-bed --out "$out/no_ambiguous"
plink --bfile "$out/no_ambiguous" --geno 0.10 --hwe 1e-10 --make-bed --out "$out/initial_variant_qc"
plink --bfile "$out/initial_variant_qc" --missing --out "$out/pre_sample_qc"
plink --bfile "$out/initial_variant_qc" --mind 0.05 --make-bed --out "$out/sample_missingness_qc"
plink --bfile "$out/sample_missingness_qc" --het --out "$out/heterozygosity"
"$R_SCRIPT" "$HIGH_HET_SCRIPT" "$out/heterozygosity.het"
# Original helper wrote a high-heterozygosity sample exclusion list.
het_remove="$out/heterozygosity.het.highhet_outside_3sd.sample.txt"
[[ -f "$het_remove" ]] || { echo "Review and supply $het_remove" >&2; exit 1; }
plink --bfile "$out/sample_missingness_qc" --remove "$het_remove" --make-bed --out "$out/heterozygosity_qc"
# Relatedness: MHC excluded for LD pruning; use prune.in, not the historical
# --exclude prune.out expression which does not restrict to the pruned subset.
awk '$1==6 && $4>24000000 && $4<36000000 {print $2}' "$out/heterozygosity_qc.bim" > "$out/mhc_exclude.snps"
plink --bfile "$out/heterozygosity_qc" --exclude "$out/mhc_exclude.snps" --maf 0.05 --indep-pairwise 50 5 0.2 --out "$out/relatedness_prune"
plink --bfile "$out/heterozygosity_qc" --extract "$out/relatedness_prune.prune.in" --genome --out "$out/relatedness"
# Historical IBD >0.9 used for near-duplicate identification; final sample
# exclusion also involved IBS/F-outlier review. Require externally reviewed list.
awk 'NR>1 && $10>0.9 {print}' "$out/relatedness.genome" > "$out/ibd_over_0.9.review.tsv"
echo 'REVIEW REQUIRED: validate strand flips, IBD pairs, IBS/F outliers, ancestry exclusions.' >&2
: "${REVIEWED_SAMPLE_REMOVE:?Set REVIEWED_SAMPLE_REMOVE to reviewed two-column FID IID exclusion file (can be empty)}"
[[ -f "$REVIEWED_SAMPLE_REMOVE" ]] || exit 1
plink --bfile "$out/heterozygosity_qc" --remove "$REVIEWED_SAMPLE_REMOVE" --make-bed --out "$out/sample_qc"
plink --bfile "$out/sample_qc" --geno 0.02 --make-bed --out "$out/variant_missingness_qc"
plink --bfile "$out/variant_missingness_qc" --mac 5 --make-bed --out "$out/mac5"
plink --bfile "$out/mac5" --hardy --freq --out "$out/reference_frequency_review"
# Historical workflow excluded chr6 markers with |frequency(target)-frequency(1KG)|>0.20
# after allele matching. This requires external ancestry-matched 1KG data.
: "${REVIEWED_CHR6_FREQ_EXCLUDE:?Set REVIEWED_CHR6_FREQ_EXCLUDE to reviewed variant-ID exclusion file (can be empty)}"
[[ -f "$REVIEWED_CHR6_FREQ_EXCLUDE" ]] || exit 1
plink --bfile "$out/mac5" --exclude "$REVIEWED_CHR6_FREQ_EXCLUDE" --make-bed --out "$out/reference_qc"
plink --bfile "$out/reference_qc" --chr 6 --from-mb 28 --to-mb 34 --make-bed --out "$out/hla_chr6_28_34Mb"
plink --bfile "$out/hla_chr6_28_34Mb" --freq --missing --out "$out/hla_chr6_28_34Mb_qc"
echo "HLA extraction complete; reference naming and strand alignment must be verified before upload."
