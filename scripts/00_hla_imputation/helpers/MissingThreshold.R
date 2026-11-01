library(ggplot2)
"%&%" <- function(a,b) paste(a,b, sep = "")

args <- commandArgs(TRUE)

# Load Data
name = args[1]
thred = args[2]
dat = read.table(name, header = TRUE)
png(name%&%".png")
ggplot(dat) +
    geom_histogram(aes(x = F_MISS, y = ..count..)) +
    theme_classic() +
    geom_vline(xintercept = as.numeric(thred), color = "red", linetype = "dashed") +
    scale_x_continuous(
        breaks = seq(0, max(dat$F_MISS, na.rm = TRUE), by = 0.005)
    ) +
    labs(x = "Genotype Missingness Rate", y = "Frequency", title = "Post 1ST SNP QC")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1,size=8))
dev.off()
