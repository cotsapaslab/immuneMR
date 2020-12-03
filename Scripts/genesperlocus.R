args <- commandArgs(TRUE)

library("biomaRt")
library(plyr)



ind <- read.table(args[1],sep="\t",header=T)
phen <- read.table(args[2],sep="\t",header=F)

phen <- phen[,1:4]
phen[,1] <- gsub("chr","",phen[,1])

#grch38 = useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")


phen$ID <- NA
for (i in 1:nrow(phen)){
phen[i,"ID"] <- unlist(strsplit(toString(phen[i,4]),"\\."))[1]} 

ann <- getBM(c("hgnc_symbol","description","chromosome_name","band","strand","start_position","end_position","ensembl_gene_id"), "ensembl_gene_id", phen$ID, grch37)
colnames(phen) <- c("CHR","START","END","GENE_ID","ensembl_gene_id")

phen <- join(phen,ann[,c("ensembl_gene_id","start_position","end_position","hgnc_symbol")],by="ensembl_gene_id",type="left")
phen$hgnc_symbol <- NULL
phen <- unique(phen)
phen <- phen[is.na(phen$start_position) == FALSE,]


for (i in 1:nrow(ind)){

chr <- ind[i,1]
snp <- as.numeric(as.character(ind[i,3]))
startbp <- as.numeric(as.character(ind[i,4]))
endbp <- as.numeric(as.character(ind[i,5]))

start <- snp-500000
end <- snp+500000

test <- phen[phen[,1] == chr,]
test$keep <- 0

test <- test[test$CHR == chr,]
test[test$start_position > start & test$start_position < end,"keep"] <- 1
test[test$start_position <= start & test$end_position > start,"keep"] <- 1
 
test <- test[test$keep == 1,]
test$keep <- NULL

output <- paste0(args[3],"PhenoList",chr,".",startbp,".",endbp,".txt")
write.table(test$GENE_ID,output,sep="\n",row.names=F,quote=F,col.names=F)



} 


