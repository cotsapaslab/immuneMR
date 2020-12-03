args <- commandArgs(TRUE)

library(plyr)


cov <- read.table(args[1],sep="\t",header=T)
PEER <- read.table(args[2],sep="",header=T, check.names=F,row.names=1)



PEER <- t(PEER)

PEER <- data.frame(PEER)
PEER$IID <- rownames(PEER)



out <- join(cov,PEER[,c(1:30,ncol(PEER))],by="IID",type="left")


out <- data.frame(out)


write.table(out,args[3],sep="\t",row.names=F,quote=F)


