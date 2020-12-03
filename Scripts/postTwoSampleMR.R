args <- commandArgs(TRUE)

library(ggplot2)
library(gridExtra)
library(ggsci)
library(plyr)

#args <- c("/home/cg859/scratch60/PRS-Project/MIP_BP/MR/TwoSampleMR/all.MRresults.txt","/home/cg859/scratch60/PRS-Project/MIP_BP/JLIMresults_selected_complete.txt","Inverse variance weighted","0.00001","/home/cg859/scratch60/PRS-Project/MIP_BP/MR/TwoSampleMR/plots_IVW_0.00001_test.png")

dat <- read.table(args[1],sep="\t",h=T)
jlim <- read.table(args[2],sep="\t",h=T)
method <- args[3]
pthres <- as.numeric(as.character(args[4]))
output <- args[5]

dat <- as.matrix(dat)
for (i in 1:nrow(dat)){dat[i,"exposure"] <- unlist(strsplit(as.character(dat[i,"exposure"]),"\\_"))[length(unlist(strsplit(as.character(dat[i,"exposure"]),"\\_")))]}
dat[dat[,"JLIM"] == "JILM+","JLIM"] <- "JLIM+"
dat <- data.frame(dat)


jlim <- as.matrix(jlim)
for (i in 1:nrow(jlim)){
if (grepl("binary",jlim[i,"maintrID"]) == FALSE){
jlim[i,"maintrID"] <- paste0(jlim[i,"maintrID"],"_norm")
}
}
jlim <- data.frame(jlim)

# Merge with JLIM data
dat$Comb <- paste0(dat$GenExpDS,"_",dat$exposure,"_",dat$outcome)
jlim$Comb <- paste0(jlim$DataSet,"_",jlim$Gene,"_",jlim$maintrID)
merged <- join(dat,jlim[,c("Comb","minp")],by="Comb",type="left")

# Select results based on method and pthreshold
res <- merged[merged$method == "MR Egger" & as.numeric(as.character(merged$pthresMR)) == pthres,]
res$pval <- as.numeric(as.character(res$pval))
res$b <- as.numeric(as.character(res$b))
res$logP <- -log10(as.numeric(as.character(res$pval)))
res <- res[is.na(res$minp) == FALSE,]
print(nrow(res))
# JLIM information
pos <- res[res$JLIM == "JLIM+",]
neg <- res[res$JLIM == "JLIM-",]
res$JLIMname <- res$JLIM
res <- as.matrix(res)
res[res[,"JLIMname"] == "JLIM+","JLIMname"] <- paste0("JLIM+","\n","(n=",as.character(nrow(pos)),")")
res[res[,"JLIMname"] == "JLIM-","JLIMname"] <- paste0("JLIM-","\n","(n=",nrow(neg),")")
res <- data.frame(res)
res$JLIMname <- as.factor(res$JLIMname)
res$pval <- as.numeric(as.character(res$pval))
res$b <- as.numeric(as.character(res$b))
res$logP <- as.numeric(as.character(res$logP))
res$absbeta <- abs(res$b)

# Comparison of pvalues between JLIM+ and JLIM- combination
pval <- wilcox.test(res[res$JLIM == "JLIM+","pval"],res[res$JLIM == "JLIM-","pval"])$p.val
pval <- formatC(pval, format = "e", digits = 2)
tit <- paste0("M.-W.-Wilcoxon test p = ",as.character(pval))
print("pvalue:")
print(tit)

p <- ggplot(res, aes(x=JLIMname, y=pval, color=JLIMname, fill=JLIMname)) + geom_boxplot(alpha=0.8,size=0.3)
p <- p + scale_color_nejm() + scale_fill_nejm()
p <- p + xlab("") + ylab("pvalue MR")
p <- p + theme(text = element_text(size=6,family="sans"), plot.title = element_text(size=5,family="sans",face="bold"), axis.text.x = element_text(family="sans",size=6,face="bold"), axis.text.y = element_text(family="sans",size=6), axis.title.x = element_text(family="sans",size=6,face = "bold"), axis.title.y = element_text(family="sans",size=6,face = "bold"))
p <- p + ggtitle(tit)
p <- p + theme_minimal()
p <- p + theme(legend.title = element_blank())
p <- p + theme(legend.position = "none")
p2 <- p

# Comparison of betas between JLIM+ and JLIM- combination
pval <- wilcox.test(res[res$JLIM == "JLIM+","b"],res[res$JLIM == "JLIM-","b"])$p.val
pval <- formatC(pval, format = "e", digits = 2)
tit <- paste0("M.-W.-Wilcoxon test p = ",as.character(pval))
print("beta:")
print(tit)

p <- ggplot(res, aes(x=JLIMname, y=absbeta, color=JLIMname, fill=JLIMname)) + geom_boxplot(alpha=0.8,size=0.3)
p <- p + scale_color_nejm() + scale_fill_nejm()
p <- p + xlab("") + ylab("Absolute estimated causal effect")
p <- p + theme(text = element_text(size=6,family="sans"), plot.title = element_text(size=5,family="sans",face="bold"), axis.text.x = element_text(family="sans",size=6,face="bold"), axis.text.y = element_text(family="sans",size=6), axis.title.x = element_text(family="sans",size=6,face = "bold"), axis.title.y = element_text(family="sans",size=6,face = "bold"))
p <- p + ggtitle(tit)
p <- p + theme_minimal()
p <- p + theme(legend.title = element_blank())
p <- p + theme(legend.position = "none")
p1 <- p

# LM function
lm_eqn <- function(df,name1,name2){
  m <- lm(y ~ x, df);
  r2=as.character(round(summary(m)$r.squared,3))
  pval=summary(m)$coefficients[2,4]
  pval <- as.character(formatC(pval, format = "e", digits = 2))
  text <- paste0("r2=",r2,"\n","p=",pval)
  return(text)
}

# Correlation JILM minp beta MR
minpJLIM=0.00001
res$minp <- as.numeric(as.character(res$minp))
res[res$minp == 0.0,"minp"] <- as.numeric(as.character(minpJLIM))
res$logminp <- -log10(res$minp)

p <- ggplot(res, aes(x=logminp,y=absbeta))  + geom_point(color="black",alpha=0.65)
p <- p + xlab("- log (minimal JLIM p value)") + ylab("Absolute estimated effect of gene expression on immune phenotype")
p <- p + geom_smooth(method=lm,color="black",level=0.95,size=0.5)
p <- p + theme(panel.background = element_rect(fill = "grey97"))
df <- res[,c("logminp","absbeta")]
colnames(df) <- c("x","y")
text <- as.character(paste("r2 =",toString(lm_eqn(df,"log minp JLIM","beta"))))
ypos <- max(res$absbeta)*0.8
xpos <- max(res$logminp)*0.8
p <- p + annotate("text",x = xpos, y= ypos,label=text)
p <- p + theme(text = element_text(size=6,family="sans"), plot.title = element_text(size=5,family="sans",face="bold"), axis.text.x = element_text(family="sans",size=6,face="bold"), axis.text.y = element_text(family="sans",size=6), axis.title.x = element_text(family="sans",size=6,face = "bold"), axis.title.y = element_text(family="sans",size=6,face = "bold"))
p <- p + theme_minimal()
p3 <- p

# Correlation JILM minp pval MR
p <- ggplot(res, aes(x=logminp,y=logP))  + geom_point(color="black",alpha=0.65)
p <- p + xlab("- log (minimal JLIM p value)") + ylab("-log (p value TwoSampleMR)")
p <- p + geom_smooth(method=lm,color="black",level=0.95,size=0.5)
p <- p + theme(panel.background = element_rect(fill = "grey97"))
df <- res[,c("logminp","logP")]
colnames(df) <- c("x","y")
text <- as.character(paste(toString(lm_eqn(df,"log minp JLIM","log MR pval"))))
ypos <- max(res$logP)*0.8
xpos <- max(res$logminp)*0.8
p <- p + annotate("text",x = xpos, y= ypos,label=text)
p <- p + theme(text = element_text(size=6,family="sans"), plot.title = element_text(size=5,family="sans",face="bold"), axis.text.x = element_text(family="sans",size=6,face="bold"), axis.text.y = element_text(family="sans",size=6), axis.title.x = element_text(family="sans",size=6,face = "bold"), axis.title.y = element_text(family="sans",size=6,face = "bold"))
p <- p + theme_minimal()
p4 <- p

# FDR p-value comparison
# FDR pvalue
res$P.adj <- p.adjust(res$pval,method="fdr",n=length(res$pval))
res$sig <- 0
res[res$pval < 0.05,"sig"] <- 1
dat <- res   
print(paste0("Number of combinations (nrow(dat)): ",nrow(dat)))
print(paste0("Number of JLIM+ combinations with FDR PRS p < 0.05: ",as.character(nrow(dat[dat$sig == 1 & dat$JLIM == "JLIM+",])),"/",nrow(dat[dat$JLIM == "JLIM+",]),"(",round(nrow(dat[dat$sig == 1 & dat$JLIM == "JLIM+",])/nrow(dat[dat$JLIM == "JLIM+",])*100,2),"%)"))
print(paste0("Number of JLIM- combinations with FDR PRS p < 0.05: ",as.character(nrow(dat[dat$sig == 1 & dat$JLIM == "JLIM-",])),"/",nrow(dat[dat$JLIM == "JLIM-",]),"(",round(nrow(dat[dat$sig == 1 & dat$JLIM == "JLIM-",])/nrow(dat[dat$JLIM == "JLIM-",])*100,2),"%)"))
print(paste0("Number of JLIM+ combinations with FDR PRS p >= 0.05: ",as.character(nrow(dat[dat$sig == 0 & dat$JLIM == "JLIM+",])),"/",nrow(dat[dat$JLIM == "JLIM+",]),"(",round(nrow(dat[dat$sig == 0 & dat$JLIM == "JLIM+",])/nrow(dat[dat$JLIM == "JLIM+",])*100,2),"%)"))
print(paste0("Number of JLIM- combinations with FDR PRS p >= 0.05: ",as.character(nrow(dat[dat$sig == 0 & dat$JLIM == "JLIM-",])),"/",nrow(dat[dat$JLIM == "JLIM-",]),"(",round(nrow(dat[dat$sig == 0 & dat$JLIM == "JLIM-",])/nrow(dat[dat$JLIM == "JLIM-",])*100,2),"%)"))







png(output,height=17,width=20,units="cm",res=300,type="cairo")
grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()
