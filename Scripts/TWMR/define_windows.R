args <- commandArgs(TRUE)


library("biomaRt")
library(plyr)



phen <- read.table(args[1],sep="\t",header=F)
exp <- read.table(args[2],sep="\t",h=T)

exp <- exp[,c(1,4)]
exp[,1] <- gsub("chr","",exp[,1])
colnames(exp) <- c("chromosome_name","phenoname")

grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

phen$ID <- NA
for (i in 1:nrow(phen)){
phen[i,"ID"] <- unlist(strsplit(toString(phen[i,1]),"\\."))[1]}

ann <- getBM(c("hgnc_symbol","description","chromosome_name","band","strand","start_position","end_position","ensembl_gene_id"), "ensembl_gene_id", phen$ID, grch37)
colnames(phen) <- c("phenoname","ensembl_gene_id")

phen <- join(phen,ann[,c("ensembl_gene_id","chromosome_name","start_position","end_position","hgnc_symbol")],by="ensembl_gene_id",type="left")
phen$hgnc_symbol <- NULL
phen <- unique(phen)
phen <- phen[is.na(phen$start_position) == FALSE,]
phen <- phen[is.na(phen$end_position) == FALSE,]

out <- phen[,c("phenoname","ensembl_gene_id","chromosome_name","start_position","end_position")]

out$window_start <- out$start_position-1000000
out$window_end <- out$end_position+1000000
out[out$window_start < 0,"window_start"] <- 0

out$chromosome_name <- NULL
out <- join(out,exp,by="phenoname",type="left")

out <- out[,c("phenoname","ensembl_gene_id","chromosome_name","start_position","end_position","window_start","window_end")]

write.table(out,args[3],sep="\t",row.names=F,quote=F,col.names=T)









