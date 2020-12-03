args <- commandArgs(TRUE)


dat1 <- read.table(args[1],sep="")
dat2 <- read.table(args[2],sep="")
output <- args[3]


dat1$Comb <- paste0(dat1$V2,"_",dat1$V5,"_",dat1$V6)
dat2$Comb <- paste0(dat2$V2,"_",dat2$V5,"_",dat2$V6)
dat2$Comb2 <- paste0(dat2$V2,"_",dat2$V6,"_",dat2$V5)

snps <- dat1[dat1$Comb %in% dat2$Comb2 == TRUE,"V2"]
write.table(snps,output,sep="\n",row.names=F,quote=F,col.names=F)






