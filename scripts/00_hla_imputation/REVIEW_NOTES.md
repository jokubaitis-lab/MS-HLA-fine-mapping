# Helper audit and change log

## Included and integrated

| File | Use | Public adaptation |
|---|---|---|
| `get_duprem_var.py` | Duplicate variant detection | Refactored; keeps lowest missingness at same position and unordered alleles; original behavior retained in principle |
| `get_remID.py` | Suggest sample to remove from IBD >0.9 pairs | Original retained; review suggestions before exclusion |
| `HighHeterozygosityOutlier.R` | High heterozygosity exclusion list | Denominator corrected to nonmissing genotype count; remove-list header disabled |
| `MissingThreshold.R` | Missingness plot | Original retained |
| `F_outlier.R` | Inbreeding F plot and removal suggestions | Original retained |
| `MissingVarThreshold.R` | Variant missingness plot | Original retained |
| `Frq1KG.R`, `Frq1KG_chr.R` | Allele frequency plots | Original retained; reference frequency preparation still requires manual integration |
| `rename_bim.py` | Match BIM variant names to 1KG | Original retained for reference, not called automatically |
