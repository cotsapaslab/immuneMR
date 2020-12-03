args <- commandArgs(TRUE)

library(plyr)


bed <- read.table(args[1],sep="\t",header=T,check.names=F)
sam <- read.table(args[2],sep="",skip=4)
fam <- read.table(args[3],sep="")
dataset <- args[5]

colnames(fam)[1:2] <- c("FID","IID")

for (j in 5:ncol(bed)){
colnames(bed)[j] <- toString(sam[j-4,1])}


rownames(bed) <- bed[,4]
bed <- bed[,5:ncol(bed)]
bed <- t(bed)

bed <- data.frame(bed)
bed$IID <- NA


if (dataset != "CMC"){
for (i in 1:nrow(bed)){
name <- rownames(bed)[i]
names <- unlist(strsplit(name,"_"))
bed[i,"IID"] <- names[length(names)]}
}else{
bed$IID <- rownames(bed)
}
out <- join(bed,fam[,1:2],by="IID",type="left")

out <- out[,c(ncol(out),(ncol(out)-1),1:(ncol(out)-2))]

write.table(out,args[4],sep="\t",row.names=F,quote=F)



