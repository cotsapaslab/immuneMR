args <- commandArgs(TRUE)

library(ggplot2)
library(ggsci)
library(plyr)
library(gridExtra)


#args <- c("/home/cg859/scratch60/PRS-Project/MIP_BP/MR/TWMR/all.alpha","/home/cg859/scratch60/PRS-Project/MIP_BP/JLIMresults_selected_complete.txt","/home/cg859/scratch60/PRS-Project/MIP_BP/test.png")

res <- read.table(args[1],sep="",h=T)
jlim <- read.table(args[2],sep="\t",h=T)

jlim <- as.matrix(jlim)
for (i in 1:nrow(jlim)){
if (grepl("binary",jlim[i,"maintrID"]) == FALSE){
jlim[i,"maintrID"] <- paste0(jlim[i,"maintrID"],"_norm")}}
jlim <- data.frame(jlim)

res$DS <- NA
res[res$DataSet == "BP","DS"] <- "Tcells"
res[res$DataSet == "BP_Mono","DS"] <- "Monocytes"
res[res$DataSet == "BP_Neutro","DS"] <- "Neutrophils"
res$Comb <- paste0(res$DS,"_",res$Gene,"_",res$ImmPhen)
jlim$Comb <- paste0(jlim$DataSet,"_",jlim$Gene,"_",jlim$maintrID)

dat <- join(res,jlim[,c("Comb","minp")],by="Comb",type="left")

# FDR pvalue
dat$P.adj <- p.adjust(dat$P,method="fdr",n=length(dat$P))
dat$sig <- 0
dat[dat$P.adj < 0.05,"sig"] <- 1
pos <- dat[dat$JLIM == "JLIM+",]
neg <- dat[dat$JLIM == "JLIM-",]
print(paste0("Number of combinations (nrow(dat)): ",nrow(dat)))
print(paste0("Number of JLIM+ combinations with FDR PRS p < 0.05: ",as.character(nrow(dat[dat$sig == 1 & dat$JLIM == "JLIM+",])),"/",nrow(dat[dat$JLIM == "JLIM+",]),"(",round(nrow(dat[dat$sig == 1 & dat$JLIM == "JLIM+",])/nrow(dat[dat$JLIM == "JLIM+",])*100,2),"%)"))
print(paste0("Number of JLIM- combinations with FDR PRS p < 0.05: ",as.character(nrow(dat[dat$sig == 1 & dat$JLIM == "JLIM-",])),"/",nrow(dat[dat$JLIM == "JLIM-",]),"(",round(nrow(dat[dat$sig == 1 & dat$JLIM == "JLIM-",])/nrow(dat[dat$JLIM == "JLIM-",])*100,2),"%)"))
print(paste0("Number of JLIM+ combinations with FDR PRS p >= 0.05: ",as.character(nrow(dat[dat$sig == 0 & dat$JLIM == "JLIM+",])),"/",nrow(dat[dat$JLIM == "JLIM+",]),"(",round(nrow(dat[dat$sig == 0 & dat$JLIM == "JLIM+",])/nrow(dat[dat$JLIM == "JLIM+",])*100,2),"%)"))
print(paste0("Number of JLIM- combinations with FDR PRS p >= 0.05: ",as.character(nrow(dat[dat$sig == 0 & dat$JLIM == "JLIM-",])),"/",nrow(dat[dat$JLIM == "JLIM-",]),"(",round(nrow(dat[dat$sig == 0 & dat$JLIM == "JLIM-",])/nrow(dat[dat$JLIM == "JLIM-",])*100,2),"%)"))


# JLIM information
dat$JLIMname <- dat$JLIM
dat <- as.matrix(dat)
dat[dat[,"JLIMname"] == "JLIM+","JLIMname"] <- paste0("JLIM+","\n","(n=",as.character(nrow(pos)),")")
dat[dat[,"JLIMname"] == "JLIM-","JLIMname"] <- paste0("JLIM-","\n","(n=",nrow(neg),")")
dat <- data.frame(dat)
dat$JLIMname <- as.factor(dat$JLIMname)


# Comparison alpha
dat$alpha <- as.numeric(as.character(dat$alpha))
dat$absalpha <- abs(dat$alpha)

pval <- wilcox.test(dat[dat$JLIM == "JLIM+","absalpha"],dat[dat$JLIM == "JLIM-","absalpha"])$p.val
pval <- formatC(pval, format = "e", digits = 2)
tit <- paste0("Mann-Whitney-Wilcoxon test p = ",as.character(pval))
print("Alpha:")
print(tit)

p <- ggplot(dat, aes(x=JLIMname, y=absalpha, color=JLIMname, fill=JLIMname)) + geom_boxplot(trim=FALSE, alpha=0.8,size=0.3)
p <- p + scale_color_nejm() + scale_fill_nejm()
p <- p + xlab("") + ylab("Absolute estimated effect of gene expression on immune phenotype")
p <- p + theme(text = element_text(size=6,family="sans"), plot.title = element_text(size=5,family="sans",face="bold"), axis.text.x = element_text(family="sans",size=6,face="bold"), axis.text.y = element_text(family="sans",size=6), axis.title.x = element_text(family="sans",size=6,face = "bold"), axis.title.y = element_text(family="sans",size=6,face = "bold"))
p <- p + ggtitle(tit)
p <- p + theme_minimal()
p <- p + theme(legend.title = element_blank()) 
p <- p + theme(legend.position = "none") 
p1 <- p

# Comparison p value
dat$logP <- -log10(as.numeric(as.character(dat$P)))

pval <- wilcox.test(dat[dat$JLIM == "JLIM+","logP"],dat[dat$JLIM == "JLIM-","logP"])$p.val
pval <- formatC(pval, format = "e", digits = 2)
tit <- paste0("Mann-Whitney-Wilcoxon test p = ",as.character(pval))
print("pvalue:")
print(tit)

p <- ggplot(dat, aes(x=JLIMname, y=logP, color=JLIMname, fill=JLIMname)) + geom_boxplot(trim=FALSE, alpha=0.8,size=0.3)
p <- p + scale_color_nejm() + scale_fill_nejm()
p <- p + xlab("") + ylab("-log(p TWMR)")
p <- p + theme(text = element_text(size=6,family="sans"), plot.title = element_text(size=5,family="sans",face="bold"), axis.text.x = element_text(family="sans",size=6,face="bold"), axis.text.y = element_text(family="sans",size=6), axis.title.x = element_text(family="sans",size=6,face = "bold"), axis.title.y = element_text(family="sans",size=6,face = "bold"))
p <- p + ggtitle(tit)
p <- p + theme_minimal()
p <- p + theme(legend.title = element_blank())
p <- p + theme(legend.position = "none") 
p2 <- p

# LM function
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

# Correlation JLIMp alpha
minpJLIM=0.00001
dat$minp <- as.numeric(as.character(dat$minp))
dat[dat$minp == 0.0,"minp"] <- as.numeric(as.character(minpJLIM))
dat$logminp <- -log10(dat$minp)

p <- ggplot(dat, aes(x=logminp,y=absalpha))  + geom_point(color="black",alpha=0.65)
p <- p + xlab("- log (minimal JLIM p value)") + ylab("Absolute estimated effect of gene expression on immune phenotype")
p <- p + geom_smooth(method=lm,color="black",level=0.95,size=0.5)
p <- p + theme(panel.background = element_rect(fill = "grey97"))
df <- dat[,c("logminp","absalpha")]
colnames(df) <- c("x","y")
text <- as.character(paste(toString(lm_eqn(df,"log minp JLIM","absolute alpha"))))
ypos <- max(dat$absalpha)*0.8
xpos <- max(dat$logminp)*0.8
p <- p + annotate("text",x = xpos, y= ypos,label=text)
p <- p + theme(text = element_text(size=6,family="sans"), plot.title = element_text(size=5,family="sans",face="bold"), axis.text.x = element_text(family="sans",size=6,face="bold"), axis.text.y = element_text(family="sans",size=6), axis.title.x = element_text(family="sans",size=6,face = "bold"), axis.title.y = element_text(family="sans",size=6,face = "bold"))
p <- p + theme_minimal()
p3 <- p

# Correlation JLIMp TWMRp

p <- ggplot(dat, aes(x=logminp,y=logP))  + geom_point(color="black",alpha=0.65)
p <- p + xlab("- log (minimal JLIM p value)") + ylab("- log (TWMR p value)")
p <- p + geom_smooth(method=lm,color="black",level=0.95,size=0.5)
p <- p + theme(panel.background = element_rect(fill = "grey97"))
df <- dat[,c("logminp","logP")]
colnames(df) <- c("x","y")
text <- as.character(paste(toString(lm_eqn(df,"log minp JLIM","log TWMR p"))))
ypos <- max(dat$logP)*0.8
xpos <- max(dat$logminp)*0.8
p <- p + annotate("text",x = xpos, y= ypos,label=text)
p <- p + theme(text = element_text(size=6,family="sans"), plot.title = element_text(size=5,family="sans",face="bold"), axis.text.x = element_text(family="sans",size=6,face="bold"), axis.text.y = element_text(family="sans",size=6), axis.title.x = element_text(family="sans",size=6,face = "bold"), axis.title.y = element_text(family="sans",size=6,face = "bold"))
p <- p + theme_minimal()
p4 <- p



# Fisher test and plot for significant TWMR results
test <- matrix(nrow=2,ncol=2)
test[1,] <- table(pos$sig)
test[2,] <- table(neg$sig)
print ("Fisher test number of sig TWMR results:")
print(test)
print(fisher.test(test))
pval <- fisher.test(test)$p.val
pval <- formatC(pval, format = "e", digits = 2)
tit <- paste0("Fisher's exact test p = ",pval)

test <- matrix(nrow=4,ncol=3,NA)
colnames(test) <- c("JLIM","TWMR.sig","Percentage")
test <- data.frame(test)
jlim0 <- paste0("JLIM-","\n","(n= ",nrow(dat[dat$JLIM == "JLIM-",]),")")
jlim1 <- paste0("JLIM+","\n","(n= ",nrow(dat[dat$JLIM == "JLIM+",]),")")
test$JLIM <- c(jlim0,jlim0,jlim1,jlim1)
test$TWMR.sig <- c("0","1","0","1")
test[test$JLIM == jlim0 & test$TWMR.sig == "0","Percentage"] <- round(100*(nrow(dat[dat$JLIM == "JLIM-" & dat$sig == 0,])/nrow(dat[dat$JLIM == "JLIM-",])),2)
test[test$JLIM == jlim1 & test$TWMR.sig == "0" ,"Percentage"] <- round(100*(nrow(dat[dat$JLIM == "JLIM+" & dat$sig == 0,])/nrow(dat[dat$JLIM == "JLIM+",])),2)
test[test$JLIM == jlim0 & test$TWMR.sig == "1","Percentage"] <- round(100*(nrow(dat[dat$JLIM == "JLIM-" & dat$sig == 1,])/nrow(dat[dat$JLIM == "JLIM-",])),2)
test[test$JLIM == jlim1 & test$TWMR.sig == "1" ,"Percentage"] <- round(100*(nrow(dat[dat$JLIM == "JLIM+" & dat$sig == 1,])/nrow(dat[dat$JLIM == "JLIM+",])),2)
test$JLIM <- as.factor(test$JLIM)
test$TWMR.sig <- c(">= 0.05","< 0.05",">= 0.05","< 0.05")
test$TWMR.sig <- as.factor(test$TWMR.sig)
test$Percentage <- as.numeric(as.character((test$Percentage)))

p <- ggplot(test, aes(x = JLIM, y = Percentage)) + geom_col(aes(fill = TWMR.sig), width = 0.3,alpha=0.8)
p <- p + scale_color_jama() + scale_fill_jama()
p <- p + xlab("") + ylab("%")
p <- p + theme(text = element_text(size=6,family="sans"), plot.title = element_text(size=5,family="sans",face="bold"), axis.text.x = element_text(family="sans",size=6,face="bold"), axis.text.y = element_text(family="sans",size=6), axis.title.x = element_text(family="sans",size=6,face = "bold"), axis.title.y = element_text(family="sans",size=6,face = "bold"))
p <- p + theme_minimal()
p <- p + ggtitle(tit)
p <- p + labs(fill = "FDR TWMR p")
p5 <- p


#Plot to output
png(args[3],width=27,height=25,units="cm",res=300,type="cairo")
grid.arrange(p2, p5, p3, p4, nrow = 2)
dev.off()


# Summary info
dat$Ngene <- as.numeric(as.character(dat$Ngene))
dat$Nsnps <- as.numeric(as.character(dat$Nsnps))
pos <- dat[dat$JLIM == "JLIM+",]
neg <- dat[dat$JLIM == "JLIM-",]
print("JLIM- Ngene (median,IQR):")
print(paste0(round(median(neg$Ngene),3),", ",round(quantile(neg$Ngene)[2],3),"-",round(quantile(neg$Ngene)[4],3)))
print("JLIM+ Ngene (median,IQR):")
print(paste0(round(median(pos$Ngene),3),", ",round(quantile(pos$Ngene)[2],3),"-",round(quantile(pos$Ngene)[4],3)))
print("JLIM- Nsnps (median,IQR):")
print(paste0(round(median(neg$Nsnps),3),", ",round(quantile(neg$Nsnps)[2],3),"-",round(quantile(neg$Nsnps)[4],3)))
print("JLIM+ Nsnps (median,IQR):")
print(paste0(round(median(pos$Nsnps),3),", ",round(quantile(pos$Nsnps)[2],3),"-",round(quantile(pos$Nsnps)[4],3)))
print("JLIM- alpha (median,IQR):")
print(paste0(round(median(neg$absalpha),3),", ",round(quantile(neg$absalpha)[2],3),"-",round(quantile(neg$absalpha)[4],3)))
print("JLIM+ alpha (median,IQR):")
print(paste0(round(median(pos$absalpha),3),", ",round(quantile(pos$absalpha)[2],3),"-",round(quantile(pos$absalpha)[4],3)))


sig <- dat[dat$sig == 1,]
nsig <- dat[dat$sig == 0,]
print("Significants Ngene (median,IQR):")
print(paste0(round(median(sig$Ngene),3),", ",round(quantile(sig$Ngene)[2],3),"-",round(quantile(sig$Ngene)[4],3)))
print("NS Ngene (median,IQR):")
print(paste0(round(median(nsig$Ngene),3),", ",round(quantile(nsig$Ngene)[2],3),"-",round(quantile(nsig$Ngene)[4],3)))
print("Significants Nsnps (median,IQR):")
print(paste0(round(median(sig$Nsnps),3),", ",round(quantile(sig$Nsnps)[2],3),"-",round(quantile(sig$Nsnps)[4],3)))
print("NS Nsnps (median,IQR):")
print(paste0(round(median(nsig$Nsnps),3),", ",round(quantile(nsig$Nsnps)[2],3),"-",round(quantile(nsig$Nsnps)[4],3)))


