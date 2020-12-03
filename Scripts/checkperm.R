args <- commandArgs(TRUE)

perm <- read.table(args[1],sep="")
nperm <- as.numeric(as.character(args[2]))
add <- read.table(args[3],sep="",header=T)
add <- add[add$TEST == "ADD",]
NAass <- add[is.na(add$P) == TRUE,]

names <- as.character(add$SNP)
colnames(perm)[2:ncol(perm)] <- names

perm <- perm[,colnames(perm) == "V1" | colnames(perm) %in% NAass$SNP == FALSE,]

perm$rem <- 0
for (i in 1:nrow(perm)){
if (TRUE %in% is.na(perm[i,1:ncol(perm)]) == TRUE){perm[i,"rem"] <- 1}}

perm[1,"rem"] <- 0
perm <- perm[perm$rem == 0,]
perm <- perm[1:(nperm+1),]
perm[,1] <- 0:nperm
perm$rem <- NULL 
write.table(perm,args[1],sep=" ",row.names=F,quote=F,col.names=F)



 
