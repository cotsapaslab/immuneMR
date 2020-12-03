args <- commandArgs(TRUE)

dat <- read.table(args[1],sep="\t")

dat$comb <- paste0(dat$V1,dat$V4,dat$V5,dat$V6)
out <- dat[duplicated(dat$comb) == FALSE,]
out$comb <- NULL

write.table(out,args[2],sep="\t",row.names=F,quote=F,col.names=F)


