#IMPUTE_recomb_rates.R
#Linearly imputes recombination rates for SNPs with missing data

args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
x <- read.table(file, header=T)

library(zoo)
x[,2] <- na.approx(x[,2],x[,1], na.rm=F)
x[is.na(x[,2]),2] <- 0

write.table(x,file,row.names = F, quote = F, eol = "\n")
