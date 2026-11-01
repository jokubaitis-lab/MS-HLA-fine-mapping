# HLA imputation and genotype QC (GRCh37/hg19)

This directory contains **expanded, public-facing reconstructions** of the genotype QC and HLA-region imputation preparation used in the study. The commands follow the original study scripts, but include explicit manual review gates where historical helper scripts, reference files, or decisions are not available. **These files are not a verified end-to-end rerun of the historical pipeline.** No individual-level genotype data, identifiers, credentials or institution-specific absolute paths are distributed.

## Two distinct imputation workflows

1. **Genome-wide SNP imputation for population structure:** the historical genome-wide workflow prepared all 22 autosomes using the HRC alignment checker, submitted data to Michigan Imputation Server, and retained genome-wide imputed SNPs for genetic PC calculation. Its historical post-imputation outputs are separate from HLA fine-mapping inputs. The precise PCA computation should be documented in the relevant PCA script.
2. **HLA-specific imputation for fine-mapping:** a separate, chr6:28–34 Mb QC and imputation workflow used the Michigan Imputation Server **Four-digit Multi-ethnic HLA v2** reference panel to obtain HLA-region dosages. Classical HLA allele and amino-acid dosage encoding for the association/fine-mapping analyses requires additional project-specific scripts.

## Files

| File | Function |
| --- | --- |
| `01_genotype_qc_mhc.sh` | Deduplication, SNPflip screening, variant/sample missingness, high-heterozygosity review, MHC-excluded LD pruning, IBD review, MAC filtering, reference-frequency review gate, HLA extraction |
| `03_prepare_hla_imputation_input.sh` | Convert a *reviewed, reference-aligned* HLA PLINK dataset to a compressed, indexed VCF |
| `02_post_imputation_qc.sh` | INFO/R2 > 0.7, dosage import, missingness summaries, sample/variant QC, MAC >= 5, VCF export |

## HLA imputation stages and historical thresholds

1. Input: GRCh37/hg19 genotype PLINK files; inspect `.bim`, `.fam`, and missingness.
2. Identify and remove duplicate variants using the historical `get_duprem_var.py` helper.
3. Assess ambiguous and strand-flipped markers with SNPflip and GRCh37 FASTA; review flipped variants separately. Palindromic removal was an optional commented step, **not** a verified executed filter.
4. Initial SNP QC: `--geno 0.10 --hwe 1e-10`.
5. Sample QC: `--mind 0.05`, high-heterozygosity outlier review (historical 3 SD helper), relatedness analysis with MHC excluded from pruning; original command history used IBD >0.9 to flag near duplicates. Additional IBS/F outlier and ancestry decisions were made in the original analysis, but their final removal lists are not included here.
6. Final SNP QC: `--geno 0.02`, remove monomorphic variants and apply `--mac 5`.
7. Compare allele frequencies against ancestry-relevant 1000 Genomes Phase 3 reference data after allele harmonisation; the historical workflow assessed absolute discrepancies >0.20 on chr6 and excluded flagged SNPs. This stage requires review of reference files and matched alleles.
8. Extract chr6:28–34 Mb (GRCh37), inspect frequencies, harmonise IDs and reference alleles; the original workflow tested `rename_bim.py`, HRC/1000G checker and its generated `Run-plink.sh`. **Do not assume every exploratory command was executed for the final submission.**
9. Upload an aligned chr6 VCF to Michigan Imputation Server with **Four-digit Multi-ethnic HLA v2** reference panel. Server-side phasing and settings should be reported from actual job metadata. Local SHAPEIT2 commands in the historical script were commented out and are not represented as executed steps.
10. Post-imputation: the original workflow explored INFO/R2 >0.5 and >0.7; the primary public script implements >0.7. Import DS dosages using PLINK2, review missingness, apply MAC >=5 and select the analysis sample set. The original code contains both 1% and 5% variant missingness settings in different branches; 5% is exposed here and should be checked against the actual analysis outputs.
11. Extract and encode HLA alleles/amino-acid positions for association and fine-mapping. **Not included:** the provided scripts do not demonstrate how HLA allele and amino-acid dosage matrices were generated; document this from the corresponding analysis scripts before claiming full reproducibility.

## Required software and resources

PLINK 1.9, PLINK2, bcftools, tabix, Bash, Python3, R, SNPflip; GRCh37 reference FASTA; 1000 Genomes Phase 3 reference data; historical `get_duprem_var.py` and high-heterozygosity helper; reviewed sample-removal and reference-frequency exclusion lists. The helper paths are provided through environment variables. The code will stop rather than silently skip unavailable QC decisions.

## Cautions and deviations from historical commands

- Relatedness pruning is expressed as `--indep-pairwise 50 5 0.2` followed by `--extract prune.in`; the historical `--indep 50 5 2` and `--exclude prune.out` combination requires review and is **not faithfully executable as written**.
- SNPflip output naming and duplicate-helper input/output conventions must be checked against installed versions before running.
- No automated replacement of `.psam` IDs based on underscores: use a verified mapping to prevent sample-ID corruption.
- INFO/R2 filtering is applied without blindly replacing VCF FILTER strings; validate against the imputation server's exact output schema.
- Reference alignment is a required manual gate, not a substitute for proper allele matching.
- The supplied files have been checked for shell syntax but **not run on individual-level data**.

## References

- [HLA QC and imputation tutorial](https://github.com/immunogenomics/HLA_analyses_tutorial/blob/main/tutorial_HLAQCImputation.ipynb)
- Michigan Imputation Server and Four-digit Multi-ethnic HLA v2 panel (specify server job settings in manuscript methods).
