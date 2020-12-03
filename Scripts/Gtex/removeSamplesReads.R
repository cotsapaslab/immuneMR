args <- commandArgs(TRUE)


met <- read.table(args[1],sep="\t",header=T)
counts <- read.table(args[2],sep="\t",skip=2,header=T,check.names=F)
rpkm <- read.table(args[3],sep="\t",skip=2,header=T,check.names=F)
tpm <- read.table(args[4],sep="\t",skip=2,header=T,check.names=F)

remove <- met[met$Mapped < 10000000,"Sample"]


counts <- counts[,colnames(counts) %in% remove == FALSE]
rpkm <- rpkm[,colnames(rpkm) %in% remove == FALSE]
tpm <- tpm[,colnames(tpm) %in% remove == FALSE]


write.table(counts,paste0(args[2],"_QC_tmp"),sep="\t",row.names=F,quote=F)
write.table(rpkm,paste0(args[3],"_QC_tmp"),sep="\t",row.names=F,quote=F)
write.table(tpm,paste0(args[4],"_QC_tmp"),sep="\t",row.names=F,quote=F)

