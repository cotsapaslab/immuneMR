args <- commandArgs(TRUE)

genes <- read.table(args[1],sep="",h=F)
dat <- read.table(args[2],sep="\t",h=T,nrows=3)

genes <- genes[genes$V1 %in% colnames(dat) == TRUE,]

write.table(genes,args[3],sep="\t",row.names=F,quote=F,col.names=F)


