args <- commandArgs(TRUE)


bim <- read.table(args[1],sep="",h=F)
ld <- read.table(args[2],sep="",h=T)
ass <- read.table(args[3],sep="",h=F)

ld <- ld[ld$R2 >= 0.2,]
bim <- bim[bim[,2] %in% ld$SNP_B == FALSE,]
ass <- ass[ass[,2] %in% bim[,2] == TRUE,]

write.table(bim,args[1],sep="\t",row.names=F,quote=F,col.names=F)
write.table(ass,args[3],sep="\t",row.names=F,quote=F,col.names=F)


