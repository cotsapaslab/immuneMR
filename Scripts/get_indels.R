args  <- commandArgs(TRUE)

bim <- read.table(args[1],sep="",header=F)
bim <- as.matrix(bim)
out <- bim[nchar(bim[,5]) > 1 | nchar(bim[,6]) > 1,]

write.table(out[,2],args[2],sep="\t",row.names=F,quote=F,col.names=F)


 
