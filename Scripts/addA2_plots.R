args <- commandArgs(TRUE)

library(plyr)

bim <- read.table(args[1],sep="")
dat <- read.table(args[2],sep="",h=T)

colnames(bim) <- c("CHR","SNP","X","BP","A1","A2")

dat <- dat[dat$SNP %in% bim$SNP == TRUE,]

out <- join(dat,bim[,c("SNP","A2")],by="SNP",type="left")

write.table(out,args[3],sep="\t",row.names=F,quote=F)

