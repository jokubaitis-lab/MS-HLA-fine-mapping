# Helper audit and change log

## Included and integrated

| File | Use | Public adaptation |
|---|---|---|
| `get_duprem_var.py` | Duplicate variant detection | Refactored; keeps lowest missingness at same position and unordered alleles; original behavior retained in principle |
| `get_remID.py` | Suggest sample to remove from IBD >0.9 pairs | Original retained; review suggestions before exclusion |
| `HighHeterozygosityOutlier.R` | High heterozygosity exclusion list | Denominator corrected to nonmissing genotype count; remove-list header disabled; **changes historical behavior** |
| `MissingThreshold.R` | Missingness plot | Original retained |
| `F_outlier.R` | Inbreeding F plot and removal suggestions | Original retained |
| `MissingVarThreshold.R` | Variant missingness plot | Original retained |
| `Frq1KG.R`, `Frq1KG_chr.R` | Allele frequency plots | Original retained; reference frequency preparation still requires manual integration |
| `rename_bim.py` | Match BIM variant names to 1KG | Original retained for reference, **not called automatically** |

## Not included in active pipeline

- `HeterozygosityOutlier.R`: superseded by high-heterozygosity-only script in historical workflow.
- `MapChrPos2Rs.py`: call was commented out in historical HLA script; includes a personal absolute path.
- `SNP2HLA.py`: separate local Beagle-based imputation implementation, not the MIS workflow. Do not conflate them.
- `beagle.jar`, `beagle2linkage.jar`, `beagle2vcf.jar`: third-party binary distributions; do not redistribute without licensing review; not needed for MIS submission.
- `HRC-1000G-check-bim.pl`, `Run-plink.sh`: not uploaded. Must be sourced/generated and reviewed separately if historical alignment is to be fully reproduced.
- `plotRsqMAF.R`, `FilterRsq.py`: referenced by the older genome-wide pipeline, not uploaded.

## Verification and security

- No uploaded private reference data, genotype data, or binary JARs included.
- Never commit output `.fam`, `.psam`, sample removal lists, imputation download links, private credentials, or individual-level data.
- Audit shell logs, exclusion counts, reference allele alignment, allele dosages and association inputs against historical runs before claiming exact manuscript reproduction.
- No validation against genotype datasets was performed. Bash syntax and Python parsing were checked only.
