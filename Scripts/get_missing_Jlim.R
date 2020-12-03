args <- commandArgs(TRUE)

res <- read.table(args[1],sep="\t",h=T)
ref <- read.table(args[2],sep="\t")

res <- res[,c("CellType","ImmPhen","Cond1stname","Locus","Gene","Cond2ndname")]

colnames(ref) <- colnames(res)
ref$Locus <- gsub("locus","locus.",ref$Locus)
res$comb <- paste0(res[,1],res[,2],res[,3],res[,4],res[,5],res[,6])
ref$comb <- paste0(ref[,1],ref[,2],ref[,3],ref[,4],ref[,5],ref[,6])


out <- ref[ref$comb %in% res$comb == FALSE,]
out$comb <- NULL
out$Locus <- gsub("locus.","locus",out$Locus)

write.table(out,args[3],sep="\t",row.names=F,quote=F,col.names=F)

