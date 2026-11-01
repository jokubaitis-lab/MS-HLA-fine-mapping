# HLA fine-mapping in multiple sclerosis

Analysis code accompanying **“HLA Fine-Mapping Resolves Multiple Sclerosis Susceptibility to Amino-Acid Residues that Shape Peptide-Binding Specificity.”**

## Overview

The major histocompatibility complex (MHC) harbours the strongest genetic associations with multiple sclerosis (MS). Extensive linkage disequilibrium and the high polymorphism of classical human leukocyte antigen (HLA) alleles make it challenging to identify the amino-acid variation underlying these associations. Fine-mapping at residue resolution provides a complementary view of HLA susceptibility by linking genetic association patterns to the peptide-binding regions of HLA molecules.

We investigated imputed classical HLA alleles and amino-acid polymorphisms in an Australian case–control dataset. Haplotype-based conditional analyses identified **HLA-DRβ1 positions 71, 13 and 86** as the principal determinants of the DRB1 association. We further evaluated residue combinations using exhaustive model searches and permutation testing, assessed their contribution to MS susceptibility using Nagelkerke's pseudo-R², and examined additional association signals at HLA class I and class II loci.

This repository provides scripts for the **genetic imputation and fine-mapping workflow** described in the manuscript. It does not distribute participant-level genotype or phenotype data. 

## Analysis workflow

The workflow follows two connected stages. Genotype QC, haplotype phasing and HLA imputation generate the classical HLA and amino-acid dosage data used in subsequent association and fine-mapping analyses. The upstream procedures were adapted to the study cohorts and analysis settings; the links below acknowledge the established tutorial on which the workflow was based.

### 1. Genotype QC and HLA imputation

**Reference:** [HLA imputation tutorial](https://github.com/immunogenomics/HLA_analyses_tutorial/blob/main/tutorial_HLAQCImputation.ipynb).

Genotypes underwent quality control, including variant and sample missingness, heterozygosity, relatedness, minor allele count and allele-frequency concordance checks. The MHC region (chr6:28–34 Mb, GRCh37/hg19) was extracted, and variants were harmonised against reference data before HLA imputation through the **Michigan Imputation Server**, using **Minimac4** and the **multi-ethnic HLA reference panel v2**. Following imputation and study-specific quality control, classical HLA alleles, amino-acid polymorphisms and MHC SNPs were prepared for association analysis.

- [`00_hla_imputation/README.md`](scripts/00_hla_imputation/README.md) summarises the genotype QC, reference harmonisation, HLA imputation and post-imputation processing stages.
- [`00_hla_imputation/01_genotype_qc_mhc.sh`](scripts/00_hla_imputation/01_genotype_qc_mhc.sh) documents pre-imputation genotype QC, sample and variant filtering, and MHC extraction.
- [`00_hla_imputation/03_prepare_hla_imputation_input.sh`](scripts/00_hla_imputation/03_prepare_hla_imputation_input.sh) prepares reference-harmonised VCF files for imputation.
- [`00_hla_imputation/02_post_imputation_qc.sh`](scripts/00_hla_imputation/02_post_imputation_qc.sh) documents post-imputation R², missingness and minor allele count filtering.
- [`00_hla_imputation/helpers/`](scripts/00_hla_imputation/helpers/) contains supporting Python and R scripts for genotype QC, outlier detection and variant harmonisation.

The scripts are adapted from the original study workflow for public sharing. Cohort-specific exclusion lists, reference datasets and server-side imputation procedures are not redistributed. Reference harmonisation and downstream HLA allele/amino-acid processing require additional external resources or study-specific steps.

### 2. HLA association and fine-mapping

**Reference:** [HLA association and fine-mapping tutorial](https://github.com/immunogenomics/HLA_analyses_tutorial/blob/main/tutorial_association.md).

The fine-mapping analyses use imputed HLA allele dosages and allele-to-residue mappings to compare logistic regression models of MS case–control status.

| Analysis | Code | Manuscript connection |
| --- | --- | --- |
| **MHC-wide single-marker association** | [`01_mhc_association/mhc_single_marker_logistic.sh`](scripts/01_mhc_association/mhc_single_marker_logistic.sh) | PLINK2 logistic association underlying the MHC association analysis|
| **DRB1 conditional fine-mapping** | [`01_drb1_conditional/conditional_test_HLA_binary_stepwise2_deviance.R`](scripts/01_drb1_conditional/conditional_test_HLA_binary_stepwise2_deviance.R) | Forward conditional haplotype-based likelihood-ratio tests identifying DRβ1 positions 71, 13 and 86|
| **Sequential MHC conditioning** | [`02_mhc_conditional/sequential_mhc_conditional.sh`](scripts/02_mhc_conditional/sequential_mhc_conditional.sh) | Adapted two-stage PLINK2 conditional testing (DRB1 classical alleles followed by MHC-wide candidates) |
| **Exhaustive residue combinations** | [`02_exhaustive/exhaustive_test.R`](scripts/02_exhaustive/exhaustive_test.R) | Comparison of residue-position pairs and triplets; 1,275 pairs and 20,825 triplets reported |
| **Permutation analysis** | [`03_permutation/permutation_test_new.R`](scripts/03_permutation/permutation_test_new.R) | Allele-to-residue mapping permutations assessing the three-position model|
| **Variance explained** | [`04_variance_explained/Variance_explained.R`](scripts/04_variance_explained/Variance_explained.R) | Nagelkerke pseudo-R² for fine-mapped residue models|
| **β86 beyond DRB1*15:01** | [`04_variance_explained/DRB1_B86_Variance_explained.R`](scripts/04_variance_explained/DRB1_B86_Variance_explained.R) | Conditional contribution of position 86 beyond the classical DRB1*15:01 association|
| **Additional MHC residue associations** | [`05_additional_mhc/univar_residue_logistic.R`](scripts/05_additional_mhc/univar_residue_logistic.R) | Residue-level association analysis at HLA loci|
| **Conditional residue/haplotype effect estimates** | [`05_additional_mhc/multiHLA_block_conditional_table.R`](scripts/05_additional_mhc/multiHLA_block_conditional_table.R) | Multivariable logistic regression of target residue haplotypes, adjusted for selected residues at other HLA loci; produces ORs, 95% CIs and frequency tables |

The table indicates the associated manuscript analyses. Publication figures require additional plotting and formatting. 

## Input data

**Individual-level study data are not publicly distributed.**

The fine-mapping scripts require appropriately formatted imputed HLA allele dosages and metadata, including:

- PLINK-style allele dosage data (`.raw`) and the corresponding allele-name file;
- an allele-to-amino-acid mapping dictionary (`HLA_DICTIONARY_AA.hg19.imgt3320.AA_tf.in_ref.rds`);
- case–control phenotypes and sample identifiers;
- ancestry principal components **PC1–PC5** and study/cohort indicators (`D_MSGP1_AU`, `D_MSGP2_AU`, `D_MSPGX_AU`), as used in the study;
- analysis-specific output prefixes and, where relevant, the HLA gene and residue-position parameters.

Allele identifiers and dosage-column ordering must be consistent across inputs. 

## Software requirements

| Stage | Software |
| --- | --- |
| Genotype processing | PLINK 1.9; standard shell utilities |
| Genotype strand/QC preparation | SNPflip and relevant reference FASTA, where applicable |
| Haplotype phasing | SHAPEIT2 |
| HLA imputation | Michigan Imputation Server / Minimac4; multi-ethnic HLA reference panel v2 |
| Fine-mapping and model comparison | R, including base `stats` functions |

The shared R scripts primarily use base R; install additional R packages if indicated by an individual script. Software versions and command-line options should be selected consistently with the input reference panel and genome build (**GRCh37/hg19** in this study).


## Citation

If you use this repository, please cite the associated manuscript:

> He J, et al. *HLA Fine-Mapping Resolves Multiple Sclerosis Susceptibility to Amino-Acid Residues that Shape Peptide-Binding Specificity.* Manuscript citation and DOI to be added upon publication.

For the foundational tutorial, please also cite the associated work indicated in [HLA_analyses_tutorial](https://github.com/immunogenomics/HLA_analyses_tutorial).

## Contacts:

**Dr. Jingni He**  
Email: [jingni.he@monash.edu](mailto:jingni.he@monash.edu)

**Prof Vilija G. Jokubaitis**  
Email: [Vilija.jokubaitis@monash.edu](mailto:Vilija.jokubaitis@monash.edu)

## Copyright License (MIT Open Source)

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
