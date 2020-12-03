args <- commandArgs(TRUE)
library(plyr)

mdsFile <- args[1]
infFile <- args[2]
pcaFile <- args[3]
sel <- args[4]
propFile <- args[5]
outputfolder <- args[6]



inf <- read.table(infFile,sep="\t",header=T) 
mds <- read.table(mdsFile,sep="",header=T)
pca <- read.table(pcaFile,sep="",header=T)
prop <- read.table(propFile,sep="\t",header=T)


if (sel == "mds"){
out <- mds[,c(1,2,4:8)]}
if (sel == "pca"){
out <- pca[,1:7]
colnames(out) <- c("FID","IID","C1","C2","C3","C4","C5")
}
out <- join(out,inf[,c("IID","AGE.V0","CMV.V0","NBYTABAC")],by="IID",type="left")
out <- join(out,prop,by="IID",type="left")

colnames(out)[8:10] <- c("Age","CMV","Tabac")

write.table(out,paste0(outputfolder,"/covariates.txt"),sep="\t",row.names=F,quote=F)


