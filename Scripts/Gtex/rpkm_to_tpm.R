args <- commandArgs(TRUE)


dat <- read.table(args[1],skip=2,sep="\t",header=T,check.names=F)

#TPM = RPKM/sum(RPKM)*1e6

for (i in 3:ncol(dat)){
dat[,i] <- as.numeric(as.character(dat[,i]))
dat[,i] <- dat[,i]/sum(dat[,i])*1000000
}

write.table(dat,args[2],sep="\t",row.names=F,quote=F)







