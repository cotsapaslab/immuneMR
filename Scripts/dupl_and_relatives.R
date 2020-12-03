args <- commandArgs(TRUE)

genome <- read.table(paste0(args[1],"_genome.genome"),h=T)
FAM <- read.table(paste0(args[1],".fam"),h=T)
IMISS <- read.table(paste0(args[1],".imiss"),h=T)

duplicates <- genome[which(genome$PI_HAT>0.7),c(1,2,3,4,10)]
duplList <- unique(c(as.character(duplicates$IID1),as.character(duplicates$IID2)))

relatives <- genome[which(genome$PI_HAT>=0.1875),c(1,2,3,4,10)]
relatives <- relatives[order(-relatives$PI_HAT),]
relList <- c(as.character(relatives$IID1),as.character(relatives$IID2))
freq <- as.data.frame(table(relList))
freq <- freq[with(freq,order(-Freq)),]

relMiss <- merge(relatives,IMISS[,c(1:2,6)],by.x=c("FID1","IID1"),by.y=c("FID","IID"))
relMiss <- merge(relMiss,IMISS[,c(1:2,6)],by.x=c("FID2","IID2"),by.y=c("FID","IID"),suffix=c("1","2"),all.x=T)
remove <- as.character(relMiss[which(relMiss$F_MISS1>=relMiss$F_MISS2),]$IID1)
remove2 <- as.character(relMiss[which(relMiss$F_MISS2>=relMiss$F_MISS1),]$IID2)
remove <- c(remove,remove2)

if (!exists("removeFinal")) removeFinal <- NULL
if (length(removeFinal)>0)
  {
  removeFinal <- c(removeFinal,remove)
  relatives <- relatives[-which(relatives$IID1 %in% removeFinal | relatives$IID2 %in% removeFinal),]
  } else if (length(remove)>0)
  {
  removeFinal <- remove
  relatives <- relatives[-which(relatives$IID1 %in% removeFinal | relatives$IID2 %in% removeFinal),]
  } else removeFinal <- remove

# Generate list of individuals to be removed
colnames(FAM) <- c("FID", "IID", "PID", "MID", "SEX", "PHENO")
removeOut <- FAM[which(FAM$IID %in% removeFinal),1:2]
write.table(removeOut, paste0(args[1],"_remove_relatives.txt",collapse=""), c=F, r=F, qu=F)


print("Number of duplicatd samples (PI_HAT > 0.7):")
print(length(duplList))
print("Number of related samples (PI_HAT > 0.1875):")
print(length(relList))





