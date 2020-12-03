args <- commandArgs(TRUE)

ass1 <- read.table(args[1],sep="",header=T) # PostJLIM/MS.chr.startbp.endbp.assoc1.txt
ass2 <- read.table(args[2],sep="",header=T) # PostJLIM/pheno.chr.startbp.endbp.assoc.linear
output <- args[3] # PostJLIM/MS.pheno.chr.startbp.endbp.lead.txt

ass1 <- ass1[is.na(ass1$P) == FALSE,]
ass1 <- ass1[ass1[,"SNP"] %in% ass2[,"SNP"] == TRUE,]
ass1 <- ass1[grepl("rs",ass1[,"SNP"]) == TRUE,]
lead <- ass1[ass1[,"P"] == min(ass1[,"P"],na.rm=T) & ass1[,"SNP"] != "NA" & is.na(ass1[,"SNP"]) == FALSE,"SNP"]
lead <- lead[is.na(lead) == FALSE]
write.table(lead,args[3],sep="\t",row.names=F,quote=F,col.names=F)


