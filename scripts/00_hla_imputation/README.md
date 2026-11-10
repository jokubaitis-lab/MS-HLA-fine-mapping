# HLA imputation

This directory contains study-adapted scripts for genotype QC, MHC preparation and post-imputation processing in the MS HLA fine-mapping analysis. Coordinates: **GRCh37/hg19**. The MIS panel selected in the historical notes was **Four-digit Multi-ethnic HLA v2**. No individual-level data, private reference panels, download tokens or executable JARs are distributed.

## Directory contents

- `QC_filtered.sh`: separate cross-cohort ANZgene–MSBase near-duplicate screening (PI_HAT > 0.98) on LD-pruned SNPs; pre-harmonised GRCh37 inputs required.
- `01_genotype_qc_mhc.sh`: initial missingness, duplicate selection, SNPflip ambiguity detection, initial `--geno 0.10 --hwe 1e-10`, sample `--mind 0.03`, high heterozygosity review, LD-pruned relatedness, suggested near-duplicate removal, inbreeding review, final `--geno 0.02 --mac 5`, reference frequency exclusion review, chr6:28–34 Mb extraction.
- `03_prepare_hla_imputation_input.sh`: converts already reference-aligned MHC PLINK data to a compressed VCF. 
- `02_post_imputation_qc.sh`: MIS dosage VCF filtering (`INFO/R2 >= 0.7`), PLINK2 dosage import and `MAC >= 5`.
- `helpers/`: selected adapted QC and plotting helpers.

## Supporting scripts

| File | Use |
|---|---|
| `get_duprem_var.py` | Duplicate variant detection |
| `get_remID.py` | Suggest sample to remove from IBD >0.9 pairs | 
| `HighHeterozygosityOutlier.R` | High and low heterozygosity exclusion list |
| `MissingThreshold.R` | Missingness plot | 
| `F_outlier.R` | Inbreeding F plot and removal suggestions |
| `MissingVarThreshold.R` | Variant missingness plot |
| `Frq1KG.R`, `Frq1KG_chr.R` | Allele frequency plots |
| `rename_bim.py` | Match BIM variant names to 1KG |

## Required software

Bash, PLINK 1.9, PLINK2, Python 3, R with `ggplot2`, SNPflip, bcftools, tabix. For alignment with an external panel, the historical pipeline used `HRC-1000G-check-bim.pl` and a generated `Run-plink.sh`; obtain those from their respective upstream distributors and inspect the generated commands before running them. Paths to external resources should be configured locally; no reference data should be committed.

## Reference resources (not provided)

- GRCh37/hg19 reference FASTA, used by SNPflip.
- 1000 Genomes Phase 3 allele frequency/genotype reference, used for concordance checking and variant naming.
- HRC r1.1 GRCh37 sites reference and the HRC-1000G checking script, used in historical alignment commands; verify applicability to the actual HLA panel.
- Four-digit Multi-ethnic HLA v2, selected in MIS. The server performs panel-specific imputation; this repository does not redistribute the panel.

## Manuscript alignment

The public QC script uses sample and final variant missingness cutoffs of 3%, and the heterozygosity helper screens both tails at 3 SD. These are manuscript-aligned adaptations, not verified historical run settings. Relatedness screening uses PI_HAT >0.9 within cohorts; cross-cohort duplicate screening (PI_HAT >0.98) was performed separately and is not implemented in this directory. The script requires manual review of sample exclusion, strand and frequency harmonisation files.
