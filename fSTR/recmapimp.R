
x <- read.table("Bot+PACcghp.s.chr22.rec", header=T)

library(zoo)
x[x$rec == 0, 2] <- NA
x[,2] <- na.approx(x[,2],x[,1], na.rm=F)
x[is.na(x[,2]),2] <- 0

write.table(x,"Bot+PACcghp.s.chr22.rec",row.names = F, quote = F, eol = "\n")
