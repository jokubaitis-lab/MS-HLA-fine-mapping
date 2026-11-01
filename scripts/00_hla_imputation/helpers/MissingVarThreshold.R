library(ggplot2)
"%&%" <- function(a,b) paste(a,b, sep = "")

args <- commandArgs(TRUE)

# Load Data
name = args[1]
#read in missingness report
dat = read.table(name, header = TRUE)
#Plot histogram of variant missingness rate 
png(paste0(name, '.png'))
ggplot(dat) +
    geom_histogram(aes(x = F_MISS, y = ..count..), binwidth = 0.01) +
    theme_classic() +
    geom_vline(xintercept = 0.02, color = "red", linetype = "dashed") +
    labs(x = "Variant Missigness Rate", y = "Frequency", title = "1st SNP and SAMPLE QC")
dev.off()

