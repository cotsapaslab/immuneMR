args <- commandArgs(TRUE)

dat <- read.table(args[1],sep="\t")
done <- read.table(args[2],sep="\t",h=T)
high <- read.table(args[3],sep="")
few <- read.table(args[4],sep="")
rep <- read.table(args[5],sep="\t")
jobid <- args[6]
out <- args[7]

high$V2 <- gsub(paste0("/home/cg859/scratch60/Logs/PRS_runJLIM_",jobid,"_"),"",high$V1)
high$V2 <- gsub("_log","",high$V2)
few$V2 <- gsub(paste0("/home/cg859/scratch60/Logs/PRS_runJLIM_",jobid,"_"),"",few$V1)
few$V2 <- gsub("_log","",few$V2)

rep$num <- 0:(nrow(rep)-1)
rep$high <- 0; rep[rep$num %in% high$V2 == TRUE,"high"] <- 1
rep$few <- 0; rep[rep$num %in% few$V2 == TRUE,"few"] <- 1


dat <- data.frame(dat)
colnames(dat) <- c("CT","IP","C1","L","G","C2")
dat$comb <- paste0(dat$CT,dat$IP,dat$C1,dat$L,dat$G,dat$C2)
colnames(rep) <- c("CT","IP","C1","L","G","C2","num","high","few")
rep$comb <- paste0(rep$CT,rep$IP,rep$C1,rep$L,rep$G,rep$C2)
done <- done[is.na(done$p) == FALSE,]
done$comb <- paste0(done$CellType,done$ImmPhen,done$Cond1stname,done$Locus,done$Gene,done$Cond2ndname); done$comb <- gsub("locus.","locus",done$comb)
dat$done <- 0; dat$rep <- 0; dat$high <- 0; dat$few <- 0; dat$Comment <- NA


# Get infos
dat[dat$comb %in% done$comb == TRUE,"done"] <- 1
dat[dat$comb %in% rep$comb == TRUE,"rep"] <- 1
dat[dat$comb %in% rep[rep$high == 1,"comb"] == TRUE,"high"] <- 1
dat[dat$comb %in% rep[rep$few == 1,"comb"] == TRUE,"few"] <- 1
dat$numrep <- NA

for (i in 1:nrow(dat)){
if (dat[i,"rep"] == 1){
dat[i,"numrep"] <- rep[rep$comb == dat[i,"comb"],"num"]}}
dat[dat$done == 0 & dat$high == 0 & dat$few == 0,"Comment"] <- "Probably too few markers"
dat[dat$high == 1,"Comment"] <- "Pmin too high"
dat[dat$few == 1,"Comment"] <- "Too few markers"

write.table(dat,out,sep="\t",row.names=F,quote=F)

print(table(dat$Comment))


