args <- commandArgs(TRUE)


#args <- c("/home/cg859/scratch60/PRS-Project/MIP_BP/PRS/PRSresults_condJLIM_5e-08.txt","/home/cg859/scratch60/PRS-Project/MIP_BP/JLIMresults_selected_complete.txt","condJLIM","0.00001","none","/home/cg859/scratch60/PRS-Project/MIP_BP/PRS/plots_PRS_5e-08_condJLIM.png","/home/cg859/scratch60/PRS-Project/MIP_BP/PRS/PRS_res_5e-08_condJLIM.txt")

library(ggplot2)
library(ggsci)
library(plyr)
library(gridExtra)

res <- read.table(args[1],sep="",h=T)
ref <- read.table(args[2],sep="",h=T)
type <- args[3]
minpJLIM <- args[4]

colnames(ref)[1] <- "ImmPhen"

res <- unique(res)

ref$locus <- paste0("locus.",ref$Chr,".",ref$StartBP,".",ref$EndBP)
res$ImmPhen <- gsub("_norm","",res$ImmPhen)

res <- data.frame(res)

if (type == "woJLIM" || type == "condJLIM"){
dat <- merge(res[,c("DataSet","Gene","ImmPhen","PRS.R2","Num_SNP","locus","P")],ref[,c("DataSet","Gene","ImmPhen","STAT","minp","locus","JLIM")],by=c("DataSet","Gene","ImmPhen","locus"))
}else if (type == "all"){
dat <- merge(res[,c("DataSet","Gene","ImmPhen","PRS.R2","Num_SNP","locus","P")],ref[,c("DataSet","Gene","ImmPhen","STAT","minp","JLIM")],by=c("DataSet","Gene","ImmPhen"))
}else{
print ("wrong type!!")
}

dat[,"minp"] <- as.numeric(as.character(dat[,"minp"]))

plus <- nrow(dat[dat$JLIM == "JLIM+",])
minus <- nrow(dat[dat$JLIM == "JLIM-",])
dat$JLIMsimple <- dat$JLIM
dat <- as.matrix(dat)
dat[dat[,"JLIM"] == "JLIM+","JLIM"] <- paste0("JLIM+","\n","(n=",plus,")")
dat[dat[,"JLIM"] == "JLIM-","JLIM"] <- paste0("JLIM-","\n","(n=",minus,")")
dat <- data.frame(dat)
dat$JLIM <- as.factor(dat$JLIM)
dat$PRS.R2 <- as.numeric(as.character(dat$PRS.R2))
dat[,"minp"] <- as.numeric(as.character(dat[,"minp"]))

# Test for difference between groups
pval <- wilcox.test(dat[dat$JLIMsimple == "JLIM+","PRS.R2"],dat[dat$JLIMsimple == "JLIM-","PRS.R2"])$p.val
pval <- formatC(pval, format = "e", digits = 2)
tit <- paste0("Mann-Whitney-Wilcoxon test p = ",as.character(pval))
print(tit)

# Number of SNPs 
dat$Num_SNP <- as.numeric(as.character(dat$Num_SNP))
print("Number of SNPs in analysis:")
print(mean(dat$Num_SNP))
print(sd(dat$Num_SNP))


# Plots
p <- ggplot(dat, aes(x=JLIMsimple, y=PRS.R2, color=JLIMsimple, fill=JLIMsimple)) + geom_violin(trim=FALSE, alpha=0.8,size=0.3)
p <- p + scale_color_nejm() + scale_fill_nejm()
p <- p + xlab("") + ylab("PRS R2")
p <- p + theme(text = element_text(size=6,family="sans"), plot.title = element_text(size=5,family="sans",face="bold"), axis.text.x = element_text(family="sans",size=6,face="bold"), axis.text.y = element_text(family="sans",size=6), axis.title.x = element_text(family="sans",size=6,face = "bold"), axis.title.y = element_text(family="sans",size=6,face = "bold"))
p <- p + ggtitle(tit)
p <- p + geom_rug(data = dat[dat$JLIMsimple == "JLIM+",],sides="r")
p <- p + geom_rug(data = dat[dat$JLIMsimple == "JLIM-",],sides="l")
p <- p + theme_minimal()
p1 <- p

lm_eqn <- function(df,name1,name2){
    m <- lm(y ~ x, df);
    print(paste0("Regression ",name1," ~ ",name2))
    print(paste0("P: ",summary(m)$coefficient[2,4]))
    print(paste0("r2: ",as.character(round(summary(m)$r.squared,3))))
 
  r2=as.character(round(summary(m)$r.squared,3))
  pval=summary(m)$coefficients[2,4]
  pval <- as.character(formatC(pval, format = "e", digits = 2))
  text <- paste0("r2=",r2,"\n","p=",pval)
  return(text)
}

dat[dat$minp == 0.0,"minp"] <- as.numeric(as.character(minpJLIM))
dat$logminp <- -log10(dat$minp)

p <- ggplot(dat, aes(x=logminp,y=PRS.R2))  + geom_point(color="black",alpha=0.65)
p <- p + xlab("- log (minimal JLIM p value)") + ylab("PRS - R2 of immune phenotype")
p <- p + geom_smooth(method=lm,color="black",level=0.95,size=0.5)
p <- p + theme(panel.background = element_rect(fill = "grey97"))
df <- dat[,c("logminp","PRS.R2")]
colnames(df) <- c("x","y")
text <- as.character(paste(toString(lm_eqn(df,"log minp JLIM","PRS R2"))))
ypos <- max(dat$PRS.R2)*0.8
xpos <- max(dat$logminp)*0.8
p <- p + annotate("text",x = xpos, y= ypos,label=text)
p <- p + theme(text = element_text(size=6,family="sans"), plot.title = element_text(size=5,family="sans",face="bold"), axis.text.x = element_text(family="sans",size=6,face="bold"), axis.text.y = element_text(family="sans",size=6), axis.title.x = element_text(family="sans",size=6,face = "bold"), axis.title.y = element_text(family="sans",size=6,face = "bold"))
p <- p + theme_minimal()
p2 <- p

# Compare based on PRS P value
dat$P <- as.numeric(as.character(dat$P))
dat$PRS.P.adj <- p.adjust(dat$P, method = "fdr", n = length(dat$P))
dat$PRS.sig <- NA
dat[is.na(dat$PRS.P.adj) == FALSE & dat$PRS.P.adj < 0.05,"PRS.sig"] <- 1
dat[is.na(dat$PRS.P.adj) == FALSE & dat$PRS.P.adj >= 0.05,"PRS.sig"] <- 0
print(paste0("Number of combinations (nrow(dat)): ",nrow(dat)))
print(paste0("Number of JLIM+ combinations with FDR PRS p < 0.05: ",as.character(nrow(dat[dat$PRS.sig == 1 & dat$JLIMsimple == "JLIM+",])),"/",nrow(dat[dat$JLIMsimple == "JLIM+",]),"(",round(nrow(dat[dat$PRS.sig == 1 & dat$JLIMsimple == "JLIM+",])/nrow(dat[dat$JLIMsimple == "JLIM+",])*100,2),"%)"))
print(paste0("Number of JLIM- combinations with FDR PRS p < 0.05: ",as.character(nrow(dat[dat$PRS.sig == 1 & dat$JLIMsimple == "JLIM-",])),"/",nrow(dat[dat$JLIMsimple == "JLIM-",]),"(",round(nrow(dat[dat$PRS.sig == 1 & dat$JLIMsimple == "JLIM-",])/nrow(dat[dat$JLIMsimple == "JLIM-",])*100,2),"%)"))
print(paste0("Number of JLIM+ combinations with FDR PRS p >= 0.05: ",as.character(nrow(dat[dat$PRS.sig == 0 & dat$JLIMsimple == "JLIM+",])),"/",nrow(dat[dat$JLIMsimple == "JLIM+",]),"(",round(nrow(dat[dat$PRS.sig == 0 & dat$JLIMsimple == "JLIM+",])/nrow(dat[dat$JLIMsimple == "JLIM+",])*100,2),"%)"))
print(paste0("Number of JLIM- combinations with FDR PRS p >= 0.05: ",as.character(nrow(dat[dat$PRS.sig == 0 & dat$JLIMsimple == "JLIM-",])),"/",nrow(dat[dat$JLIMsimple == "JLIM-",]),"(",round(nrow(dat[dat$PRS.sig == 0 & dat$JLIMsimple == "JLIM-",])/nrow(dat[dat$JLIMsimple == "JLIM-",])*100,2),"%)"))
# Fisher test

tryCatch({
testmat <- as.matrix(table(dat[,c("PRS.sig","JLIMsimple")]))
print("Fisher test PRS significance for JLIM+ and JILM- combinations:")
fisher.test(testmat)$p.value
pval <- fisher.test(testmat)$p.value
pval <- formatC(pval, format = "e", digits = 2)
tit <- paste0("Fisher test p = ",as.character(pval))
print(tit)
}, error = function(e) {
})

# Plot bar plot
test <- matrix(nrow=4,ncol=3,NA)
colnames(test) <- c("JLIM","PRS.sig","Percentage")
test <- data.frame(test)
sigtit0 <- paste0(">= 0.05","\n","(n= ",nrow(dat[dat$PRS.sig == 0,]),")")
sigtit1 <- paste0("< 0.05","\n","(n= ",nrow(dat[dat$PRS.sig == 1,]),")")
test$PRS.sig <- c(sigtit0,sigtit0,sigtit1,sigtit1)
test$JLIM <- c("-","+","-","+")
test[test$PRS.sig == sigtit0 & test$JLIM == "-","Percentage"] <- round(100*(nrow(dat[dat$JLIM == "JLIM-" & dat$PRS.sig == 0,])/nrow(dat[dat$JLIM == "JLIM-",])),2)
test[test$PRS.sig == sigtit1 & test$JLIM == "-" ,"Percentage"] <- round(100*(nrow(dat[dat$JLIM == "JLIM-" & dat$PRS.sig == 1,])/nrow(dat[dat$JLIM == "JLIM-",])),2)
test[test$PRS.sig == sigtit0 & test$JLIM == "+","Percentage"] <- round(100*(nrow(dat[dat$JLIM == "JLIM+" & dat$PRS.sig == 0,])/nrow(dat[dat$JLIM == "JLIM+",])),2)
test[test$PRS.sig == sigtit1 & test$JLIM == "+" ,"Percentage"] <- round(100*(nrow(dat[dat$JLIM == "JLIM+" & dat$PRS.sig == 1,])/nrow(dat[dat$JLIM == "JLIM+",])),2)
jlimname0 <- as.character(dat[dat$JLIM == "JLIM-","JLIM"][1])
jlimname1 <- as.character(dat[dat$JLIM == "JLIM+","JLIM"][1])
test$JLIM <- c(jlimname0,jlimname1,jlimname0,jlimname1)
test$JLIM <- as.factor(test$JLIM)
test$PRS.sig <- as.factor(test$PRS.sig)
test$Percentage <- as.numeric(as.character((test$Percentage)))
#colnames(test)[2] <- PRS.pvalue

p <- ggplot(test, aes(x = JLIM, y = Percentage)) + geom_col(aes(fill = PRS.sig), width = 0.3,alpha=0.8)
p <- p + scale_color_jama() + scale_fill_jama()
p <- p + xlab("") + ylab("%")
p <- p + theme(text = element_text(size=6,family="sans"), plot.title = element_text(size=5,family="sans",face="bold"), axis.text.x = element_text(family="sans",size=6,face="bold"), axis.text.y = element_text(family="sans",size=6), axis.title.x = element_text(family="sans",size=6,face = "bold"), axis.title.y = element_text(family="sans",size=6,face = "bold"))
p <- p + theme_minimal()
p <- p + ggtitle(tit)
p3 <- p


# Plot correlation JLIM minp to PRS P
dat$logPRSP <- -log10(dat$P)
#png(args[6],width=8,height=7,units="cm",res=300,type="cairo")
p <- ggplot(dat, aes(x=logminp,y=logPRSP))  + geom_point(color="black",alpha=0.65)
p <- p + xlab("- log (minimal JLIM p value)") + ylab("- log (PRS p value)")
p <- p + geom_smooth(method=lm,color="black",level=0.95,size=0.5)
p <- p + theme(panel.background = element_rect(fill = "grey97"))
df <- dat[,c("logminp","logPRSP")]
colnames(df) <- c("x","y")
text <- as.character(paste(toString(lm_eqn(df,"log minp JLIM","log PRS p"))))
ypos <- max(dat$logPRSP,na.rm=T)*0.8
xpos <- max(dat$logminp,na.rm=T)*0.8
p <- p + annotate("text",x = xpos, y= ypos,label=text)
p <- p + theme(text = element_text(size=6,family="sans"), plot.title = element_text(size=5,family="sans",face="bold"), axis.text.x = element_text(family="sans",size=6,face="bold"), axis.text.y = element_text(family="sans",size=6), axis.title.x = element_text(family="sans",size=6,face = "bold"), axis.title.y = element_text(family="sans",size=6,face = "bold"))
p <- p + theme_minimal()
p4 <- p
#dev.off()

print(table(dat$JLIMsimple))

print("JLIM+:")
print(mean(dat[dat$JLIMsimple == "JLIM+","PRS.R2"],na.rm=T))
print("JLIM-:")
print(mean(dat[dat$JLIMsimple == "JLIM-","PRS.R2"],na.rm=T))


png(args[6],width=27,height=25,units="cm",res=300,type="cairo")
grid.arrange(p1, p3, p2, p4, nrow = 2)
dev.off()


# Save results for next step (check JLIM-/PRS+ pairs)
dat <- dat[order(dat$PRS.R2,decreasing=T),]
dat$JLIM <- NULL
write.table(dat,args[7],sep="\t",row.names=F,quote=F)







