library(ggplot2)
"%&%" <- function(a,b) paste(a,b, sep = "")

args <- commandArgs(TRUE)

#read in file
name = args[1]
#name = 'hgdp_chr6.1stSNPQC.1stSampQC.miss.frq.diff.plot'
dat = read.table(name)
#make an AF by AF plot comparing target and 1KG datasets
png(paste0(name, '.png'))
ggplot(dat, aes(x = V2, y = V3)) +
    geom_point() +
    theme_classic() +
    labs(x = "Target Dataset Freq", y = "1KG Freq") +
    geom_abline(intercept = c(0.20, -0.20), size = 1.5, linetype = "dashed", color = "red")
dev.off()

