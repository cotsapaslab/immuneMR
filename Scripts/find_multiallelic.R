args <- commandArgs(TRUE)

print(args[2])
print(args[1])

bim <- read.table(args[1],sep="")
dup <- bim[duplicated(bim[,2]) == TRUE,2]
dups <- bim[bim[,2] %in% dup == TRUE,]
write.table(dups[,2],args[2],row.names=F,quote=F,col.names=F)



