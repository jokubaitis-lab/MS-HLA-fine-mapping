##run R script to identify outliers
##This script plots the F inbreeding coeeficients (het from IBS)
"%&%" <- function(a,b) paste(a,b, sep = "")

args <- commandArgs(TRUE)

# Load Data
name = args[1]
out=name
## load fies
het <- read.table(name, header=T)
png(paste(out, "png", sep="."))
plot(het$F, main="Inbreeding coefficients analysis", xlab="", ylab="Inbreeding coefficient F")
dev.off()
## list of individuals to exclude
## cut-off set to +-0.05 (n=14)
ex1 <- subset(het, F>=0.05 | F<=(-0.05), select=c(FID, IID))
write.table(ex1, paste(name, "exc0_05", "remove", sep="."), row.names=F, col.names=F, quote=F)

