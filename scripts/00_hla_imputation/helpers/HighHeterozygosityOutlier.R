# Heterozygosity = (N.NM. - O.HOM.) / N.NM. (PLINK --het).
# The historical helper used a different denominator; sample flags may differ.
# Both tails (>3 SD from cohort mean) are screened as described in Methods.
library(ggplot2)
"%&%" <- function(a,b) paste(a,b, sep = "")

args <- commandArgs(TRUE)

#read in file
name = args[1]
dat = read.table(name, header = TRUE)
#calculate heterozygosity rate per indvidual
dat[, 'het'] = (dat[, 'N.NM.'] - dat[, 'O.HOM.'])/dat[, 'N.NM.']
m = mean(dat[, 'het'])
sd = sd(dat[, 'het'])
#Create sample ID column that matches fam naming convention
dat[, 'sample'] = paste(dat[, 'FID'], dat[, 'IID'], sep = '_')
#Plotting heterozygosity distribution as a histogram
png(paste0(name, '.png'))
ggplot(dat) +
    geom_histogram(aes(x = het, y = ..count..), binwidth = 0.002) +
    geom_vline(xintercept = c(m - (3 * sd), m, m + (3 * sd)), color = "red", linetype = "dashed") +
    theme_classic() +
    labs(x = "Heterozygosity per Sample", y = "Frequency", title = "Outlier Heterozygosity Identification")
dev.off()

high = dat[!is.na(dat[, 'het']) & (dat[, 'het'] > m + (3 * sd) | dat[, 'het'] < m - (3 * sd)), ]
write.table(high[, c(1:2)], paste0(name, ".highhet_outside_3sd.sample.txt"), quote=F, col.names=F, row.names=F)

