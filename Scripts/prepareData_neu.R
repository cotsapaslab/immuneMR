args <- commandArgs(TRUE)


dat <- read.table(args[1],sep="",h=T)
bim <- read.table(args[2],sep="",h=F)
freq <- read.table(args[3],sep="",h=T)

pheno <- args[4]
output <- args[5]

library(plyr)
library(gCMAP)

colnames(bim)[2] <- "SNP"

    dat$Z <- zScores(dat$P)
    dat$SE <- abs(dat$BETA/dat$Z)
    dat$Z <- NULL

dat <- plyr::join(dat,bim[,c(2,5,6)],by="SNP",type="left")

dat$OA <- NA
for (i in 1:nrow(dat)){
if (dat[i,"V5"] == dat[i,"A1"]){dat[i,"OA"] <- as.character(dat[i,"V6"])
}else{
dat[i,"OA"] <- as.character(dat[i,"V5"])} 
}

dat <- dat[,c("SNP","A1","OA","P","BETA","SE")]

freq <- freq[freq$SNP %in% dat$SNP == TRUE,]
freq$comb <- paste0(freq$SNP,"_",freq$A1)
dat$comb <- paste0(dat$SNP,"_",dat$A1)
dat1 <- plyr::join(dat,freq[,c("comb","MAF")],by="comb",type="left")
dat <- dat1
dat$comb <- NULL

colnames(dat) <- c("SNP","effect_allele","other_allele","pval","beta","se","eaf") 
dat$Phenotype <- pheno
write.table(dat,output,sep="\t",row.names=F,quote=F)





