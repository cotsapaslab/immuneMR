args <- commandArgs(TRUE)

library(plyr)
library(gCMAP)

assfile <- args[1]
bimfile <- args[2]
freqfile <- args[3]
pheno <- args[4]
output <- args[5]

bim <- read.table(bimfile,sep="",h=F)
colnames(bim)[2] <- "SNP"
freq <- read.table(freqfile,sep="",h=T)
  
  # Read assoc and calculate beta and se
  dat <- read.table(assfile,sep="",h=T)
  
  if ("OR" %in% colnames(dat) == TRUE){
    dat$Z <- zScores(dat$P)
    dat$BETA <- log(dat$OR)
    dat$SE <- abs(dat$BETA/dat$Z)
    dat$Z <- NULL
    dat <- plyr::join(dat,bim[,c(2,5,6)],by="SNP",type="left")
    dat$OA <- NA
    dat[dat$A1 == dat$V5,"OA"] <- as.character(dat[dat$A1 == dat$V5,"V6"])
    dat[dat$A1 == dat$V6,"OA"] <- as.character(dat[dat$A1 == dat$V6,"V5"])
    dat <- dat[,c("SNP","A1","OA","P","BETA","SE")]


  }else{
    dat$Z <- zScores(dat$P)
    dat$SE <- abs(dat$BETA/dat$Z)
    dat$Z <- NULL
    dat <- plyr::join(dat,bim[,c(2,5,6)],by="SNP",type="left")
    dat$OA <- NA
    dat[dat$A1 == dat$V5,"OA"] <- as.character(dat[dat$A1 == dat$V5,"V6"])
    dat[dat$A1 == dat$V6,"OA"] <- as.character(dat[dat$A1 == dat$V6,"V5"]) 
    dat <- dat[,c("SNP","A1","OA","P","BETA","SE")]

}

freq <- freq[freq$SNP %in% dat$SNP == TRUE,]
freq$comb <- paste0(freq$SNP,"_",freq$A1)
dat$comb <- paste0(dat$SNP,"_",dat$A1)
dat1 <- plyr::join(dat,freq[,c("comb","MAF")],by="comb",type="left")
dat <- dat1
dat$comb <- NULL


colnames(dat) <- c("SNP","effect_allele","other_allele","pval","beta","se","eaf")
dat$Phenotype <- pheno
write.table(dat,output,sep="\t",row.names=F,quote=F)

