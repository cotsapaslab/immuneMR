args <- commandArgs(TRUE)



dat <- read.table(args[1],sep="",h=T)
datall <- read.table(args[2],sep="",h=T)
first <- args[3]
outdir <- args[4]

print("read")
#Remove MHC region
dat <- dat[dat$CHR != 6 | (dat$BP < 25000000 | dat$BP > 35000000),]
datall <- datall[datall$CHR != 6 | (datall$BP < 25000000 | datall$BP > 35000000),]

for (i in 1:nrow(dat)){
print(i)
lead <- dat[i,]
leadSNP <- as.character(lead[1,"SNP"])
leadpos <- lead[1,"BP"]
leadCHR <- lead[1,"CHR"]

start <- leadpos - 100000
end <- leadpos + 100000

out <- datall[datall$CHR == leadCHR & datall$BP >= start & datall$BP <= end,]

write.table(out,paste0(outdir,"/",first,".",leadCHR,".",start,".",end,".txt"),sep="\t",row.names=F,quote=F)

out <- matrix(nrow=1,ncol=5)
out[1,1] <- leadCHR
out[1,2] <- leadSNP
out[1,3] <- leadpos
out[1,4] <- start
out[1,5] <- end
colnames(out) <- c("CHR","SNP","BP","STARTBP","ENDBP")

indFile <- paste0(outdir,"/",first,"-index.tsv")

if(!file.exists(indFile)){
write.table(out,paste0(outdir,"/",first,"-index.tsv"),sep="\t",row.names=F,quote=F)
}else{
text <- paste0(out[1,1],"\t",out[1,2],"\t",out[1,3],"\t",out[1,4],"\t",out[1,5])
write(text,sep="\t",file=indFile,append=TRUE)
}
}


