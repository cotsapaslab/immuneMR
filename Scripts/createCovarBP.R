args <- commandArgs(TRUE)

library(plyr)

fam <- read.table(args[1],sep="")
inf <- read.table(args[2],sep="\t",header=T)
mds <- read.table(args[3],sep="",h=T)
pca <- read.table(args[4],sep="",h=T)
sel <- args[5]

colnames(fam)[1:2] <- c("FID","IID")
colnames(inf)[1] <- c("IID")

if (sel == "mds"){
mds <- mds[,c(1,2,4:8)]
}else{
mds <- pca[,1:7]
colnames(mds) <- c("FID","IID","C1","C2","C3","C4","C5")
}

out <- join(fam[,1:2],inf[,c("IID","DONOR_AGE")],by="IID",type="left")  
out <- join(out,mds,by=c("FID","IID"),type="left")

out$Age <- NA
out[out$DONOR_AGE == "20 - 25","Age"] <- 1
out[out$DONOR_AGE == "25 - 30","Age"] <- 2
out[out$DONOR_AGE == "30 - 35","Age"] <- 3
out[out$DONOR_AGE == "35 - 40","Age"] <- 4
out[out$DONOR_AGE == "40 - 45","Age"] <- 5
out[out$DONOR_AGE == "45 - 50","Age"] <- 6
out[out$DONOR_AGE == "50 - 55","Age"] <- 7
out[out$DONOR_AGE == "55 - 60","Age"] <- 8
out[out$DONOR_AGE == "60 - 65","Age"] <- 9
out[out$DONOR_AGE == "65 - 70","Age"] <- 10
out[out$DONOR_AGE == "70 - 75","Age"] <- 11
out[out$DONOR_AGE == "75 - 80","Age"] <- 12
out[out$DONOR_AGE == "80 - 85","Age"] <- 13
out$DONOR_AGE <- NULL

write.table(out,args[6],sep="\t",row.names=F,quote=F)



