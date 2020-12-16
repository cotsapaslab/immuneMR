args <- commandArgs(TRUE)

library(tidyr)

dat <- read.table(args[1],sep="\t",h=T)

dat <- dat[,c("ImmPhen","CellType","Gene","Cond1stname","Cond1st","Cond2ndname","Cond2nd","idxSNP","idxBP","idxP","minP2","STAT","p","Locus")]
dat <- unique(dat)

# Remove unconditioned association signals when there are conditioned
dat$Comb <- paste0(dat$CellType,"_",dat$ImmPhen,"_",dat$Cond1st,"_",dat$Locus,"_",dat$Gene)
conds <- dat[is.na(dat$Cond2nd) != TRUE,]
dat$Conditioned2nd <- 0
dat[dat$Comb %in% conds$Comb == TRUE,"Conditioned2nd"] <- 1
dat$Comb <- paste0(dat$CellType,"_",dat$ImmPhen,"_",dat$Cond2nd,"_",dat$Locus,"_",dat$Gene)
conds <- dat[is.na(dat$Cond1st) != TRUE,]
dat$Conditioned1st <- 0
dat[dat$Comb %in% conds$Comb == TRUE,"Conditioned1st"] <- 1
dat <- dat[(is.na(dat$Cond2nd) != TRUE | dat$Conditioned2nd == 0) & (is.na(dat$Cond1st) != TRUE | dat$Conditioned1st == 0),]
dat$Comb <- NULL

dat$Comb <- paste0(dat$ImmPhen,"_",dat$CellType,"_",dat$Gene)


print("Total number of perfromed JLIM tests:" )
print(nrow(dat))
print("Total number of unique CT/ImmPhen/Gene combinations:")
print(length(unique(dat$Comb)))
print("Total number of unique ImmPhens:")
print(length(unique(dat$ImmPhen)))
print("Total number of unqiue genes:")
print(length(unique(dat$Gene)))
print("Total nuber of unique loci:")
print(length(unique(dat$Locus)))


# Calulcate FDR p values
dat$FDRp <- p.adjust(dat$p,method="fdr",n=nrow(dat))
dat$JLIM <- NA
dat$JLIMFDR <- NA
dat[dat$p < 0.05,"JLIM"] <- "JLIM+"
dat[dat$p >= 0.05,"JLIM"] <- "JLIM-"
dat[dat$FDRp < 0.05,"JLIMFDR"] <- "JLIM+"
dat[dat$FDRp >= 0.05,"JLIMFDR"] <- "JLIM-"


dat <- separate(dat,col="Locus",into=c("LN","chrom","StartBP","EndBP"),sep="\\.",remove=F)

print("JLIM+ vs. JLIM- total:")
print(table(dat$JLIMFDR))

neg <- dat[dat$JLIMFDR == "JLIM-",]
pos <- dat[dat$JLIMFDR == "JLIM+",]

# Select only one JLIM results per unique ImmPhen/CT/Gene combiation
out <- dat[,c("ImmPhen","CellType","Gene","Comb")]; out <- unique(out)
out$chrom <- NA; out$StartBP <- NA; out$EndBP <- NA; out$idxSNP <- NA; out$idxBP <- NA; out$STAT <- NA; out$minp <- NA; out$Cond <- NA; out$JLIM <- NA; out$JLIMFDR <- NA; out$Cond1stname <- NA; out$Cond2ndname <- NA
for (i in 1:nrow(out)){
if (out[i,"Comb"] %in% pos$Comb == TRUE){
	#JLIM + combinations
	test <- pos[pos$Comb == out[i,"Comb"] & pos$p == min(pos[pos$Comb == out[i,"Comb"],"p"],na.rm=T),]
	if (nrow(test) > 2){print(paste0(i,"!!!!!!!!!!!!!!!"))}
	test$STAT <- abs(test$STAT)
	test <- test[order(test$STAT,decreasing=T),]
	names <- c("chrom","StartBP","EndBP","idxSNP","idxBP","STAT","JLIM","JLIMFDR")
	for (n in 1:length(names)){name <- names[n]; out[i,name] <- test[1,name]}		
	out[i,"minp"] <- test[1,"p"]; out[i,"Cond"] <- as.character(test[1,"Cond2nd"]);out[i,"idxSNP"] <- as.character(test[1,"idxSNP"])
	out[i,"Cond1stname"] <- as.character(test[1,"Cond1stname"]); out[i,"Cond2ndname"] <- as.character(test[1,"Cond2ndname"])
}else{
	# JLIM - combinations
	test <- neg[neg$Comb == out[i,"Comb"] & neg$p == min(neg[neg$Comb == out[i,"Comb"],"p"],na.rm=T),]
	anzahl <- nrow(test)
	randnum <- sample(1:anzahl,1)
	names <- c("chrom","StartBP","EndBP","idxSNP","idxBP","STAT","JLIM","JLIMFDR")
        for (n in 1:length(names)){name <- names[n]; out[i,name] <- test[randnum,name]}
        out[i,"minp"] <- test[randnum,"p"]; out[i,"Cond"] <- as.character(test[randnum,"Cond2nd"]); out[i,"idxSNP"] <- as.character(test[randnum,"idxSNP"])
	out[i,"Cond1stname"] <- as.character(test[randnum,"Cond1stname"]); out[i,"Cond2ndname"] <- as.character(test[randnum,"Cond2ndname"])
}
}


print("Number of unique cell type / gene / immPhen / locus combinations:")
print(nrow(out))
print("JLIM+ and JLIM- combinations:")
print(table(out$JLIMFDR))

out$Comb <- NULL

write.table(out,args[2],sep="\t",row.names=F,quote=F)
write.table(dat,args[3],sep="\t",row.names=F,quote=F)
