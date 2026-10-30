# Genotype quality control and HLA imputation

This directory documents the genotype preparation and HLA imputation stages underlying the study of HLA amino-acid variation and multiple sclerosis susceptibility. The analysis used GRCh37/hg19 coordinates, genotype QC, HLA-region preparation, and imputation using the Michigan Imputation Server with the four-digit multi-ethnic HLA reference panel v2. See the [HLA QC and imputation tutorial](https://github.com/immunogenomics/HLA_analyses_tutorial/blob/main/tutorial_HLAQCImputation.ipynb) for reference harmonisation and detailed platform-specific procedures.

## Scripts

| Script | Purpose |
| --- | --- |
| [`01_genotype_qc_mhc.sh`](01_genotype_qc_mhc.sh) | Generates initial missingness summaries, applies representative variant filters and extracts chr6:28–34 Mb for HLA preparation. |
| [`02_post_imputation_qc.sh`](02_post_imputation_qc.sh) | Filters a Michigan HLA dosage VCF by INFO/R2 > 0.7, imports dosages with PLINK2 and applies MAC >= 5. |

These are **adapted, parameterised excerpts of the study workflow**, not a standalone implementation of every cohort-specific QC or server-side imputation step. The initial script does not perform duplicate removal, SNPflip strand checks, sample outlier/relatedness exclusion, ancestry selection, or allele-frequency harmonisation; those steps were performed upstream in the historical workflow and require cohort-specific reference resources and review. The imputation service and reference-panel settings must be configured in the Michigan Imputation Server; no access links or credentials are included here.

## Study workflow

1. Assess variant and sample missingness, duplicated variants, strand alignment, heterozygosity, relatedness, ancestry and reference allele-frequency concordance.
2. Extract the HLA interval on chromosome 6 (28–34 Mb, GRCh37) and prepare reference-aligned imputation input.
3. Submit compatible input to the Michigan Imputation Server with the four-digit multi-ethnic HLA reference panel v2. SHAPEIT2 phasing is described in the study methods; the original command history also contains a commented local SHAPEIT2 branch, so it is not presented here as a verified executed command.
4. Retain variants meeting imputation quality R2 > 0.7 and minor allele count >= 5; prepare aligned HLA allele and amino-acid dosages for fine-mapping.

## Usage

```bash
bash scripts/00_hla_imputation/01_genotype_qc_mhc.sh /path/to/genotypes /path/to/qc_output
bash scripts/00_hla_imputation/02_post_imputation_qc.sh /path/to/chr6.dose.vcf.gz /path/to/hla_output
```

Requirements: PLINK 1.9, PLINK2, bcftools, tabix and Bash; additional original QC tools include SNPflip, Python, R, 1000 Genomes reference resources and helper scripts from the upstream tutorial. Individual-level genotype and phenotype data are not distributed.
