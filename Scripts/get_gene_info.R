args <- commandArgs(TRUE)


library("biomaRt")
library(plyr)
library(tidyr)


grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

dat <- read.table(args[1],sep="",h=F)
dat <- separate(data=dat,col="V1",into=c("GeneID","remove"),sep="\\.",remove=F)

ann <- getBM(c("chromosome_name","start_position","ensembl_gene_id"), "ensembl_gene_id", dat$GeneID, grch37)

colnames(dat)[2] <- "ensembl_gene_id"

out <- merge(dat,ann,by="ensembl_gene_id")


out <- out[,c("V1","chromosome_name","start_position")]
write.table(out,args[2],sep=" ",row.names=F,quote=F,col.names=F)

