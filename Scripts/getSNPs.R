args <- commandArgs(TRUE)
library(plyr)

bim1 <- read.table(args[1],sep="")
bim2 <- read.table(args[2],sep="")
outputdir <- args[3]
colnames(bim1) <- c("CHR","SNP","X","BP","A1","A2")
colnames(bim2) <- colnames(bim1)

# Save original bim-Files
orig1 <- bim1
orig2 <- bim2

# Functions

#:::::::
#get other allele (A/T, C/G...)
#:::::::

other <- function(al){
if (al =="A"){return("T")
}else if (al =="T"){return("A")
}else if (al =="C"){return("G")
}else if (al =="G"){return("C")
}}

# Create combinations of BP; CHR, and alleles for comparison
bim1$comb <- paste0(bim1$CHR,"_",bim1$BP)
bim2$comb <- paste0(bim2$CHR,"_",bim2$BP)
bim1$comb1 <- paste0(bim1$CHR,"_",bim1$BP,"_",bim1$A1,"_",bim1$A2)
bim1$comb2 <- paste0(bim1$CHR,"_",bim1$BP,"_",bim1$A2,"_",bim1$A1)
bim2$comb1 <- paste0(bim2$CHR,"_",bim2$BP,"_",bim2$A1,"_",bim2$A2)
bim2$comb2 <- paste0(bim2$CHR,"_",bim2$BP,"_",bim2$A2,"_",bim2$A1)

# Only keep biallelic SNPs
multi1 <- bim1[nchar(as.character(bim1$A2)) > 1 | nchar(as.character(bim1$A2)) > 1,]
multi2 <- bim2[nchar(as.character(bim2$A2)) > 1 | nchar(as.character(bim2$A2)) > 1,]
bim1 <- bim1[nchar(as.character(bim1$A1)) == 1 & nchar(as.character(bim1$A2)) == 1,]
bim2 <- bim2[nchar(as.character(bim2$A1)) == 1 & nchar(as.character(bim2$A2)) == 1,]

# Only select SNPs present in both data sets (by position)
bim1 <- bim1[bim1$comb %in% bim2$comb == TRUE,]
bim2 <- bim2[bim2$comb %in% bim1$comb == TRUE,] 

# Check alleles
same1 <- bim1[bim1$comb1 %in% bim2$comb1 == TRUE,]
same2 <- bim2[bim2$comb1 %in% bim1$comb1 == TRUE,]
efff1 <- bim1[bim1$comb1 %in% bim2$comb2 == TRUE,]
efff2 <- bim2[bim2$comb1 %in% bim1$comb2 == TRUE,]
diff1 <- bim1[bim1$comb1 %in% bim2$comb1 == FALSE & bim1$comb1 %in% bim2$comb2 == FALSE,]
diff2 <- bim2[bim2$comb1 %in% bim1$comb1 == FALSE & bim2$comb1 %in% bim1$comb2 == FALSE,]
if (nrow(diff1) > 0){
for (i in 1:nrow(diff1)){
diff1[i,"A1f"] <- other(diff1[i,"A1"])
diff1[i,"A2f"] <- other(diff1[i,"A2"])}
diff1$comb1f <- paste0(diff1$CHR,"_",diff1$BP,"_",diff1$A1f,"_",diff1$A2f)
diff1$comb2f <- paste0(diff1$CHR,"_",diff1$BP,"_",diff1$A2f,"_",diff1$A1f)
diff1$flip <- 0
diff1[diff1$comb1f %in% bim2$comb1 == TRUE | diff1$comb1f %in% bim2$comb2 == TRUE,"flip"] <- 1
}
if (nrow(diff2) > 0){
for (i in 1:nrow(diff2)){
diff2[i,"A1f"] <- other(diff2[i,"A1"])
diff2[i,"A2f"] <- other(diff2[i,"A2"])}
diff2$comb1f <- paste0(diff2$CHR,"_",diff2$BP,"_",diff2$A1f,"_",diff2$A2f)
diff2$comb2f <- paste0(diff2$CHR,"_",diff2$BP,"_",diff2$A2f,"_",diff2$A1f)
diff2$flip <- 0
diff2[diff2$comb1f %in% bim1$comb1 == TRUE | diff2$comb1f %in% bim1$comb2 == TRUE,"flip"] <- 1
}
remdiff1 <- diff1[diff1$flip == 0,]
remdiff2 <- diff2[diff2$flip == 0,]
flip1 <- diff1[diff1$flip == 1,]
flip2 <- diff2[diff2$flip == 1,]

# Get amibguous SNPs
amb1 <- bim1[(bim1$A1 == "A" & bim1$A2 == "T") | (bim1$A1 == "T" & bim1$A2 == "A") | (bim1$A1 == "C" & bim1$A2 == "G") | (bim1$A1 == "G" & bim1$A2 == "C"),]  
amb2 <- bim2[(bim2$A1 == "A" & bim2$A2 == "T") | (bim2$A1 == "T" & bim2$A2 == "A") | (bim2$A1 == "C" & bim2$A2 == "G") | (bim2$A1 == "G" & bim2$A2 == "C"),]

# Get duplicates
dups1 <- bim1[duplicated(bim1$comb) == TRUE,]
dups2 <- bim2[duplicated(bim2$comb) == TRUE,]
dup1 <- bim1[bim1$comb %in% dups1$comb == TRUE,]
dup2 <- bim2[bim2$comb %in% dups2$comb == TRUE,]
bim1 <- bim1[bim1$comb %in% dups1$comb == FALSE,]
bim2 <- bim2[bim2$comb %in% dups2$comb == FALSE,] 
bim1 <- bim1[bim1$comb %in% bim2$comb == TRUE,]
bim2 <- bim2[bim2$comb %in% bim1$comb == TRUE,]

# Use names from bim2
colnames(orig1) <- c("CHR.1","SNP.1","X.1","BP.1","A1.1","A2.1")
orig1$comb <- paste0(orig1$CHR.1,"_",orig1$BP.1)
end2 <- bim2[(bim2$comb1 %in% same2$comb1 == TRUE & bim2$comb2 %in% same2$comb2 == TRUE) | (bim2$comb1 %in% efff2$comb1 == TRUE & bim2$comb2 %in% efff2$comb2 == TRUE),] 
end2 <- end2[end2$comb %in% dups2$comb == FALSE,]
updtnames <- join(orig1,end2,by="comb",type="left")
updtnames <- as.matrix(updtnames)
updtnames[is.na(updtnames[,"SNP"]) == TRUE,"SNP"] <- as.character(updtnames[is.na(updtnames[,"SNP"]) == TRUE,"SNP.1"])

# Get update names for flips
colnames(flip1)[1:7] <- c("CHR.1","SNP.1","X.1","BP.1","A1.1","A2.1","comb")
colnames(flip2)[1:7] <- c("CHR.2","SNP.2","X.2","BP.2","A1.2","A2.2","comb")
updtnames2 <- join(flip1,flip2,by="comb",type="left")

# SNPs to keep 
write.table(bim1$SNP,paste0(outputdir,"/keep1.txt"),sep="\n",row.names=F,quote=F,col.names=F)
write.table(bim2$SNP,paste0(outputdir,"/keep2.txt"),sep="\n",row.names=F,quote=F,col.names=F)

#Output lists
write.table(remdiff1$SNP,paste0(outputdir,"/remove_diff1.txt"),sep="\n",row.names=F,quote=F,col.names=F)
write.table(remdiff2$SNP,paste0(outputdir,"/remove_diff2.txt"),sep="\n",row.names=F,quote=F,col.names=F)
write.table(flip1$SNP,paste0(outputdir,"/flip1.txt"),sep="\n",row.names=F,quote=F,col.names=F)
write.table(flip2$SNP,paste0(outputdir,"/flip2.txt"),sep="\n",row.names=F,quote=F,col.names=F)
write.table(multi1$SNP,paste0(outputdir,"/multi1.txt"),sep="\n",row.names=F,quote=F,col.names=F)
write.table(multi2$SNP,paste0(outputdir,"/multi2.txt"),sep="\n",row.names=F,quote=F,col.names=F)
write.table(amb1$SNP,paste0(outputdir,"/remove_amb1.txt"),sep="\n",row.names=F,quote=F,col.names=F)
write.table(amb2$SNP,paste0(outputdir,"/remove_amb2.txt"),sep="\n",row.names=F,quote=F,col.names=F)
write.table(dup1$SNP,paste0(outputdir,"/remove_dup1.txt"),sep="\n",row.names=F,quote=F,col.names=F)
write.table(dup2$SNP,paste0(outputdir,"/remove_dup2.txt"),sep="\n",row.names=F,quote=F,col.names=F)
write.table(updtnames[,c("SNP.1","SNP")],paste0(outputdir,"/updatenames.txt"),sep="\t",row.names=F,quote=F,col.names=F)
write.table(updtnames2[,c("SNP.1","SNP.2")],paste0(outputdir,"/updatenamesflips.txt"),sep="\t",row.names=F,quote=F,col.names=F)



