# HLA amino-acid fine-mapping in multiple sclerosis

## Overview

This repository provides the analysis scripts accompanying our study of HLA amino-acid variation and multiple sclerosis (MS) susceptibility in an Australian case–control cohort (4,114 cases and 2,631 controls). The study investigates how specific residues within HLA peptide-binding regions account for association signals traditionally attributed to classical HLA alleles.

The central analysis identifies three HLA-DRB1 positions—β71, β13 and β86—and evaluates their combined contribution to MS susceptibility. Additional analyses examine the robustness of the three-position model, variation explained by HLA residues, and association signals outside DRB1.

**Scope:** This repository shares the statistical analysis code and its relationship to the reported findings. Individual-level genotype and phenotype data are not distributed. Structural modelling and proteomic analyses are outside the scope of this repository.

## Analysis workflow

| Step | Script | Manuscript analysis / reported result |
| --- | --- | --- |
| 1. DRB1 conditional fine-mapping | [`conditional_test_HLA_binary_stepwise2_deviance.R`](scripts/01_drb1_conditional/conditional_test_HLA_binary_stepwise2_deviance.R) | Forward conditional haplotype-based likelihood-ratio tests; selection of DRB1 β71, β13 and β86 |
| 2. Exhaustive combinations | [`exhaustive_test.R`](scripts/02_exhaustive/exhaustive_test.R) | Pairwise and three-position model search (1,275 pairs; 20,825 triplets reported) |
| 3. Permutation analysis | [`permutation_test_new.R`](scripts/03_permutation/permutation_test_new.R) | Empirical significance of the three-position model (10,000 permutations; reported P = 0.0046) |
| 4. Variance explained | [`Variance_explained(1).R`](scripts/04_variance_explained/Variance_explained(1).R) | Nagelkerke pseudo-R² for DRB1 residue models and HLA association models |
| 5. β86 conditional contribution | [`DRB1_B86_Variance_explained(1).R`](scripts/04_variance_explained/DRB1_B86_Variance_explained(1).R) | Contribution of β86 beyond DRB1*15:01 |
| 6. Additional MHC associations | [`univar_residue_logistic.R`](scripts/05_additional_mhc/univar_residue_logistic.R) | Residue-level association follow-up outside DRB1; not a complete substitute for all multi-locus conditional analyses |

The scripts produce association and model-comparison tables. Final manuscript figure layouts may require separate plotting code.

## Input data

The analysis scripts operate on preprocessed imputed HLA allele dosages and associated metadata. Depending on the script, required inputs include:

- PLINK-format allele dosage export (`*.raw`) and the matching allele-name file;
- an HLA allele-to-amino-acid dictionary (`HLA_DICTIONARY_AA.hg19.imgt3320.AA_tf.in_ref.rds`);
- case–control phenotype file with individual identifiers;
- covariate file with PC1–PC5 and cohort indicators (`D_MSGP1_AU`, `D_MSGP2_AU`, `D_MSPGX_AU`);
- an output prefix and, where applicable, a target gene and excluded positions.

**Data format and allele ordering must be consistent across the dosage matrix and allele-name file.** Input data are not included in this repository.

## Software

- **R** with packages `data.table`, `dplyr`, `tidyr`, and `stringr` (requirements vary by script).
- Standard R `stats` functions for binomial logistic regression, likelihood-ratio tests and model comparisons.

Install the main packages with:

```r
install.packages(c("data.table", "dplyr", "tidyr", "stringr"))
```

The code was developed for an R-based HLA fine-mapping workflow. Exact package versions depend on the original analysis environment.

## Running the scripts

Scripts use positional command-line arguments; see the `args <- commandArgs(TRUE)` section at the beginning of each script for the expected order. For example, the exhaustive search accepts the dictionary, dosage prefix, allele-name file, phenotype and covariate files, gene, maximum model size and output prefix:

```bash
Rscript scripts/02_exhaustive/exhaustive_test.R \
  /path/to/hla_dictionary.rds \
  /path/to/hla_dosages \
  /path/to/allele_names.txt \
  /path/to/phenotypes.txt MS_pheno \
  /path/to/covariates.txt \
  "PC1 + PC2 + PC3 + PC4 + PC5 + D_MSGP1_AU + D_MSGP2_AU + D_MSPGX_AU" \
  DRB1 3 /path/to/results/drb1_exhaustive
```

For the DRB1 stepwise script, set `HLA_DICTIONARY_DIR` to the directory containing the dictionary RDS before running it. The public DRB1 conditional script consistently uses PC1–PC5 and all three cohort indicators across rounds.

## Interpretation and citation

The code is provided to document the analytical methods underlying the manuscript. Without the study's individual-level input data, the repository is not a standalone reproduction of the published numerical results. The reported statistics above refer to the manuscript, not outputs bundled here.

**Citation:** Please cite the associated HLA fine-mapping manuscript when using this workflow. Add the final article DOI and citation upon publication.
