args <- commandArgs(TRUE)

library(plyr)

dat <- read.table(args[1],sep="\t",header=T)
fam <- read.table(args[2],sep="")


dat$DONOR_HEALTH_STATUS <- NA
dat$DISEASE <- NA
dat$donor_id <- NA
dat$DONOR_AGE <- NA
dat$MOLECULE <- NA
dat$BIOMATERIAL_PROVIDER <- NA
dat$BIOMATERIAL_TYPE <- NA
dat$TISSUE_TYPE <- NA
dat$SAMPLE_ONTOLOGY <- NA
dat$phenotype <- NA
dat$DONOR_SEX <- NA
dat$DONOR_ETHNICITY <- NA
dat$TISSE_DEPOT <- NA
dat$gender <- NA
dat$"ENA-CHECKLIST" <- NA
dat$SAMPLE_ONTOLOGY_URI <- NA
dat$TISSUE_DEPOT <- NA
dat$COLLECTION_METHOD <- NA
dat$DISEASE_ONTOLOGY_URI <- NA

for (i in 1:nrow(dat)){
infos <- unlist(strsplit(as.character(dat[i,1]),";"))
for (j in 1:length(infos)){
info <- unlist(strsplit(infos[j],"="))[1]
par <- unlist(strsplit(infos[j],"="))[2]
if (info %in% colnames(dat) == TRUE){
dat[i,info] <- par
}
}}

dat[,1] <- dat$donor_id
colnames(dat)[1] <- "DONOR_ID"
dat[,"donor_id"] <- NULL


dat <- unique(dat)

write.table(dat,args[3],sep="\t",row.names=F,quote=F)


colnames(fam)[2] <- "DONOR_ID"
fam <- join(fam,dat[,c("DONOR_ID","gender")],by="DONOR_ID",type="left") 
fam[fam$gender == "male",5] <- 1
fam[fam$gender == "female",5] <- 2
fam$gender <- NULL

write.table(fam,args[2],sep=" ",row.names=F,quote=F,col.names=F)


