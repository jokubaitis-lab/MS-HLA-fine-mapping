# HLA imputation — public workflow (review draft)

This repository reconstructs the **HLA-region QC and Michigan Imputation Server (MIS)** workflow from historical scripts for an MS HLA amino-acid fine-mapping study. Coordinates: **GRCh37/hg19**. The MIS panel selected in the historical notes was **Four-digit Multi-ethnic HLA v2**. No individual-level data, private reference panels, download tokens or executable JARs are distributed.

**Status: review draft, not yet an end-to-end validated reproduction.** It is important to check the historical execution logs and final sample/variant counts before calling this a manuscript-reproduction pipeline. In particular, the historical QC script had inconsistent or experimental branches, and the final HLA allele/amino-acid dosage export is not supplied.

## Directory contents

- `01_genotype_qc_mhc.sh`: initial missingness, duplicate selection, SNPflip ambiguity detection, initial `--geno 0.10 --hwe 1e-10`, sample `--mind 0.05`, high heterozygosity review, LD-pruned relatedness, suggested near-duplicate removal, inbreeding review, final `--geno 0.02 --mac 5`, reference frequency exclusion review, chr6:28–34 Mb extraction.
- `03_prepare_hla_imputation_input.sh`: converts **already reference-aligned** MHC PLINK data to a compressed VCF. **Does not implement allele alignment.**
- `02_post_imputation_qc.sh`: MIS dosage VCF filtering (`INFO/R2 > 0.7`), PLINK2 dosage import and `MAC >= 5`.
- `helpers/`: selected adapted QC and plotting helpers. See `REVIEW_NOTES.md` for changes and excluded files.

## Required software

Bash, PLINK 1.9, PLINK2, Python 3, R with `ggplot2`, SNPflip, bcftools, tabix. For alignment with an external panel, the historical pipeline used `HRC-1000G-check-bim.pl` and a generated `Run-plink.sh`; obtain those from their respective upstream distributors and **inspect the generated commands** before running them. Paths to external resources should be configured locally; no reference data should be committed.

## Example workflow (only after review)

```bash
# Paths to data and reference FASTA are local to the analyst.
export REVIEWED_SAMPLE_REMOVE=/private/verified_sample_removals.txt
export REVIEWED_CHR6_FREQ_EXCLUDE=/private/verified_chr6_frequency_exclusions.txt
bash 01_genotype_qc_mhc.sh /private/genotypes/study /private/results/qc /private/references/human_g1k_v37.fasta

# STOP: review flagged strand mismatches, genotype IDs, ancestry and frequency harmonisation.
# Apply and verify reference alignment outside this script before proceeding.
REFERENCE_ALIGNMENT_CONFIRMED=yes bash 03_prepare_hla_imputation_input.sh /private/results/aligned_mhc /private/results/upload/hla

# Submit hla.vcf.gz to MIS with the study-selected HLA panel; download chr6.dose.vcf.gz.
bash 02_post_imputation_qc.sh /private/results/chr6.dose.vcf.gz /private/results/hla_postqc
```

The sample removal and chr6 frequency-exclusion files are **manual review gates**, not automatically inferred. An empty file is valid only when review establishes no exclusions are needed. Avoid publishing sample IDs or output QC reports that might contain individual-level information.

## Reference resources (not redistributed)

- GRCh37/hg19 reference FASTA, used by SNPflip.
- 1000 Genomes Phase 3 allele frequency/genotype reference, used for concordance checking and variant naming.
- HRC r1.1 GRCh37 sites reference and the Will Rayner HRC-1000G checking script, used in historical alignment commands; verify applicability to the actual HLA panel.
- Four-digit Multi-ethnic HLA v2, selected in MIS. The server performs panel-specific imputation; this repository does not redistribute the panel.

## Distinct genome-wide workflow

A separate historical genome-wide MIS SNP imputation pipeline (chromosomes 1–22, HRC-oriented alignment) supported the calculation of genome-wide genetic PCs. It is **not** the same as HLA-specific imputation and is not reproduced by these scripts. PC computation itself is not provided here.

## Caveats

- `SNP2HLA.py` and Beagle JARs were supplied as additional files but were **not present as executed steps in the later MIS HLA workflow**. They are not included in the public pipeline pending confirmation of their use, licensing and dependencies.
- The supplied `rename_bim.py` is included as a reviewed helper but is not invoked automatically: matching chr:position and unordered alleles is insufficient for safe REF/ALT alignment.
- The supplied `get_remID.py` suggests removals for pairs with IBD PI_HAT > 0.9; final relatedness decisions require review.
- `HighHeterozygosityOutlier.R` has a corrected denominator compared with the uploaded historical file; **this may change results** and must be checked against the executed analysis.
- The original post-imputation workflow also explored R2 > 0.5 and additional missingness and sample-ID steps. The minimal public script only implements the R2 > 0.7 and MAC >= 5 path, not every experimental branch.
