args <- commandArgs(TRUE)

library(ggplot2)
library(ggsci)
library(plyr)
library(gridExtra)
library(grid)


fancy_scientific <- function(l) {
     # turn in to character string in scientific notation
     l <- format(l, scientific = TRUE)
     # quote the part before the exponent to keep all the digits
     l <- gsub("^(.*)e", "'\\1'e", l)
     # turn the 'e+' into plotmath format
     l <- gsub("e", "%*%10^", l)
     # return this as an expression
     parse(text=l)
}

#:::::::::::::::
# Part 1 
#:::::::::::::::

res <- read.table(args[1],sep="",h=T)
ref <- read.table(args[2],sep="",h=T)
type <- args[3]
minpJLIM <- args[4]

colnames(ref)[1] <- "ImmPhen"

ref$locus <- paste0("locus.",ref$Chr,".",ref$StartBP,".",ref$EndBP)
res$ImmPhen <- gsub("_norm","",res$ImmPhen)

res <- data.frame(res)
res <- unique(res)

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
#dat[dat[,"JLIM"] == "JLIM+","JLIM"] <- paste0("JLIM+","\n","(n=",plus,")")
#dat[dat[,"JLIM"] == "JLIM-","JLIM"] <- paste0("JLIM-","\n","(n=",minus,")")
dat[dat[,"JLIM"] == "JLIM+","JLIM"] <- paste0("shared effect \n (n=",plus,")")
dat[dat[,"JLIM"] == "JLIM-","JLIM"] <- paste0("no shared effect \n (n=",minus,")")
dat <- data.frame(dat)
dat$JLIM <- as.factor(dat$JLIM)
dat$PRS.R2 <- as.numeric(as.character(dat$PRS.R2))
dat[,"minp"] <- as.numeric(as.character(dat[,"minp"]))

# Test for difference between groups
pval <- wilcox.test(dat[dat$JLIMsimple == "JLIM+","PRS.R2"],dat[dat$JLIMsimple == "JLIM-","PRS.R2"])$p.val
pval <- formatC(pval, format = "e", digits = 2)
#pval <- fancy_scientific(pval)
tit <- paste0("MWW p=",as.character(pval))
print(tit)

# Number of SNPs 
dat$Num_SNP <- as.numeric(as.character(dat$Num_SNP))
print("Number of SNPs in analysis:")
print(mean(dat$Num_SNP))
print(sd(dat$Num_SNP))


# Plots
p <- ggplot(dat, aes(x=JLIM, y=PRS.R2, color=JLIM, fill=JLIM)) + geom_violin(trim=FALSE, alpha=0.8,size=0.3)
p <- p + scale_color_nejm() + scale_fill_nejm()
p <- p + xlab("") + ylab("PRS.R2")
p <- p + geom_rug(data = dat[dat$JLIMsimple == "JLIM+",],sides="r")
p <- p + geom_rug(data = dat[dat$JLIMsimple == "JLIM-",],sides="l")
p <- p + theme_minimal()
p <- p + theme(legend.position = "none")
p <- p + theme(panel.grid = element_blank())
p <- p + theme(plot.title = element_text(size=7.5,family="sans",face="bold"), axis.text.x = element_text(family="sans",size=7.5,face="bold"), axis.text.y = element_text(family="sans",size=7.5), axis.title.x = element_text(family="sans",size=7.5,face = "bold"), axis.title.y = element_text(family="sans",size=7.5,face = "bold"))
xpos <- 1.25
ypos <- max(dat$PRS.R2)*0.9
text <- tit
p <- p + ylab(bquote(bold(R[PRS]^2)))
#p <- p + annotate("text",x = xpos, y= ypos,label=text,fontface=2,size=2,hjust=0)
p <- p + ggtitle(tit)
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

p <- ggplot(dat, aes(x=logminp,y=PRS.R2))  + geom_point(color="grey50",alpha=0.65)
p <- p + xlab(bquote(bold('-log'[10]~'(JLIM p)'))) + ylab("PRS.R2")
p <- p + geom_smooth(method=lm,color="black",level=0.95,size=0.5)
p <- p + theme(panel.background = element_rect(fill = "grey97"))
df <- dat[,c("logminp","PRS.R2")]
colnames(df) <- c("x","y")
text <- as.character(paste(toString(lm_eqn(df,"log minp JLIM","PRS R2"))))
ypos <- max(dat$PRS.R2)*0.8
ypos2 <- max(dat$PRS.R2)*0.72
xpos <- max(dat$logminp)*0.65
p <- p + ylab(bquote(bold(R[PRS]^2)))
p <- p + theme_minimal()
p <- p + theme(panel.grid = element_blank())
p <- p + theme(plot.title = element_text(size=7.5,family="sans",face="bold"), axis.text.x = element_text(family="sans",size=7.5,face="bold"), axis.text.y = element_text(family="sans",size=7.5), axis.title.x = element_text(family="sans",size=7.5,face = "bold"), axis.title.y = element_text(family="sans",size=7.5,face = "bold"))
p <- p + annotate("text",x = xpos, y= ypos,label=text,fontface=2,size=2.5,hjust=0)
p2 <- p

# Compare based on PRS P value
dat$P <- as.numeric(as.character(dat$P))
dat$PRS.P.adj <- p.adjust(dat$P, method = "fdr", n = length(dat$P))
dat$PRS.sig <- NA
dat[is.na(dat$PRS.P.adj) == FALSE & dat$PRS.P.adj < 0.05,"PRS.sig"] <- 1
dat[is.na(dat$PRS.P.adj) == FALSE & dat$PRS.P.adj >= 0.05,"PRS.sig"] <- 0
#dat[is.na(dat$PRS.P.adj) == FALSE & dat$P < 0.05,"PRS.sig"] <- 1
#dat[is.na(dat$PRS.P.adj) == FALSE & dat$P >= 0.05,"PRS.sig"] <- 0

# Fisher test

tryCatch({
testmat <- as.matrix(table(dat[,c("PRS.sig","JLIM")]))
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
test[test$PRS.sig == sigtit0 & test$JLIM == "-","Percentage"] <- round(100*(nrow(dat[dat$JLIMsimple == "JLIM-" & dat$PRS.sig == 0,])/nrow(dat[dat$JLIMsimple == "JLIM-",])),2)
test[test$PRS.sig == sigtit1 & test$JLIM == "-" ,"Percentage"] <- round(100*(nrow(dat[dat$JLIMsimple == "JLIM-" & dat$PRS.sig == 1,])/nrow(dat[dat$JLIMsimple == "JLIM-",])),2)
test[test$PRS.sig == sigtit0 & test$JLIM == "+","Percentage"] <- round(100*(nrow(dat[dat$JLIMsimple == "JLIM+" & dat$PRS.sig == 0,])/nrow(dat[dat$JLIMsimple == "JLIM+",])),2)
test[test$PRS.sig == sigtit1 & test$JLIM == "+" ,"Percentage"] <- round(100*(nrow(dat[dat$JLIMsimple == "JLIM+" & dat$PRS.sig == 1,])/nrow(dat[dat$JLIMsimple == "JLIM+",])),2)
jlimname0 <- as.character(dat[dat$JLIMsimple == "JLIM-","JLIM"][1])
jlimname1 <- as.character(dat[dat$JLIMsimple == "JLIM+","JLIM"][1])
test$JLIM <- c(jlimname0,jlimname1,jlimname0,jlimname1)
test$JLIM <- as.factor(test$JLIM)
test$PRS.sig <- as.factor(test$PRS.sig)
test$Percentage <- as.numeric(as.character((test$Percentage)))
#colnames(test)[2] <- PRS.pvalue

p <- ggplot(test, aes(x = JLIM, y = Percentage)) + geom_col(aes(fill = PRS.sig), width = 0.3,alpha=0.8)
p <- p + scale_color_jama() + scale_fill_jama(name="FDR PRS p",labels=c("<0.05",">=0.05"))
p <- p + xlab("") + ylab("%")
p <- p + theme(text = element_text(size=6,family="sans"), plot.title = element_text(size=5,family="sans",face="bold"), axis.text.x = element_text(family="sans",size=6,face="bold"), axis.text.y = element_text(family="sans",size=6), axis.title.x = element_text(family="sans",size=6,face = "bold"), axis.title.y = element_text(family="sans",size=6,face = "bold"))
p <- p + theme_minimal()
#p <- p + ggtitle(tit)
p <- p + theme(plot.title = element_text(size=7.5,family="sans",face="bold"), axis.text.x = element_text(family="sans",size=7.5,face="bold"), axis.text.y = element_text(family="sans",size=7.5), axis.title.x = element_text(family="sans",size=7.5,face = "bold"), axis.title.y = element_text(family="sans",size=7.5,face = "bold"))
#p <- p + annotate("text",x = xpos, y= ypos,label=text,fontface=2,size=2,hjust=0)
p <- p + theme(legend.title=element_text(size=7.5,face="bold"),legend.text = element_text(size=7.5))
p <- p + theme(legend.key.size = unit(0.3, "cm"),legend.margin=margin(0,0,0,0),legend.box.margin=margin(0,-5,0,-30))
p <- p + theme(axis.title.y = element_text(vjust = -3))
p3 <- p


# Plot correlation JLIM minp to PRS P
dat$logPRSP <- -log10(dat$P)
#png(args[6],width=8,height=7,units="cm",res=300,type="cairo")
p <- ggplot(dat, aes(x=logminp,y=logPRSP))  + geom_point(color="grey50",alpha=0.65)
p <- p + xlab("-log(JLIM p)") + ylab("-log(PRS p)")
p <- p + geom_smooth(method=lm,color="black",level=0.95,size=0.5)
p <- p + theme(panel.background = element_rect(fill = "grey97"))
df <- dat[,c("logminp","logPRSP")]
colnames(df) <- c("x","y")
text <- as.character(paste(toString(lm_eqn(df,"log minp JLIM","log PRS p"))))
ypos <- max(dat$logPRSP,na.rm=T)*0.9
xpos <- max(dat$logminp,na.rm=T)*0.7
p <- p + theme_minimal()
p <- p + theme(plot.title = element_text(size=7.5,family="sans",face="bold"), axis.text.x = element_text(family="sans",size=6,face="bold"), axis.text.y = element_text(family="sans",size=6), axis.title.x = element_text(family="sans",size=6,face = "bold"), axis.title.y = element_text(family="sans",size=6,face = "bold"))
p <- p + annotate("text",x = xpos, y= ypos,label=text,fontface=2,size=2,hjust=0)
p4 <- p
#dev.off()




#:::::::::::::::
# Part 1 
#:::::::::::::::

res <- read.table(args[5],sep="",h=T)
ref <- read.table(args[6],sep="",h=T)
type <- args[7]
minpJLIM <- args[8]

colnames(ref)[1] <- "ImmPhen"

ref$locus <- paste0("locus.",ref$Chr,".",ref$StartBP,".",ref$EndBP)
res$ImmPhen <- gsub("_norm","",res$ImmPhen)

res <- data.frame(res)
res <- unique(res)

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
#dat[dat[,"JLIM"] == "JLIM+","JLIM"] <- paste0("JLIM+","\n","(n=",plus,")")
#dat[dat[,"JLIM"] == "JLIM-","JLIM"] <- paste0("JLIM-","\n","(n=",minus,")")
dat[dat[,"JLIM"] == "JLIM+","JLIM"] <- paste0("shared effect\n(n=",plus,")")
dat[dat[,"JLIM"] == "JLIM-","JLIM"] <- paste0("no shared effect\n(n=",minus,")")
dat <- data.frame(dat)
dat$JLIM <- as.factor(dat$JLIM)
dat$PRS.R2 <- as.numeric(as.character(dat$PRS.R2))
dat[,"minp"] <- as.numeric(as.character(dat[,"minp"]))

# Test for difference between groups
pval <- wilcox.test(dat[dat$JLIMsimple == "JLIM+","PRS.R2"],dat[dat$JLIMsimple == "JLIM-","PRS.R2"])$p.val
pval <- formatC(pval, format = "e", digits = 2)
#tit <- paste0("Mann-Whitney-Wilcoxon test p = ",as.character(pval))
tit <- paste0("MWW p=",as.character(pval))
print(tit)

# Number of SNPs 
dat$Num_SNP <- as.numeric(as.character(dat$Num_SNP))
print("Number of SNPs in analysis:")
print(mean(dat$Num_SNP))
print(sd(dat$Num_SNP))


# Plots
#dat$logPRS.R2 <- log10(dat$PRS.R2)
p <- ggplot(dat, aes(x=JLIM, y=PRS.R2, color=JLIM, fill=JLIM)) + geom_violin(trim=FALSE, alpha=0.8,size=0.3)
p <- p + scale_color_nejm() + scale_fill_nejm()
p <- p + xlab("") + ylab("PRS.R2")
p <- p + geom_rug(data = dat[dat$JLIMsimple == "JLIM+",],sides="r")
p <- p + geom_rug(data = dat[dat$JLIMsimple == "JLIM-",],sides="l")
p <- p + theme_minimal()
p <- p + theme(panel.grid = element_blank())
xpos <- 1.25
ypos <- max(dat$PRS.R2)*0.9
text <- tit
#p <- p + annotate("text",x = xpos, y= ypos,label=text,fontface=2,size=2,hjust=0)
p <- p + ggtitle(tit)
p <- p + ylab(bquote(bold(R[PRS]^2)))
p <- p + theme(legend.position = "none")
p <- p + theme(plot.title = element_text(size=7.5,family="sans",face="bold"), axis.text.x = element_text(family="sans",size=7.5,face="bold"), axis.text.y = element_text(family="sans",size=7.5), axis.title.x = element_text(family="sans",size=7.5,face = "bold"), axis.title.y = element_text(family="sans",size=7.5,face = "bold"))
p5 <- p

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

p <- ggplot(dat, aes(x=logminp,y=PRS.R2))  + geom_point(color="grey50",alpha=0.65)
p <- p + xlab(bquote(bold('-log'[10]~'(JLIM p)'))) + ylab("PRS.R2")
p <- p + geom_smooth(method=lm,color="black",level=0.95,size=0.5)
p <- p + theme(panel.background = element_rect(fill = "grey97"))
df <- dat[,c("logminp","PRS.R2")]
colnames(df) <- c("x","y")
text <- as.character(paste(toString(lm_eqn(df,"log minp JLIM","PRS R2"))))
ypos <- max(dat$PRS.R2)*0.9
ypos2 <- max(dat$PRS.R2)*0.82
xpos <- max(dat$logminp)*0.65
p <- p + theme_minimal()
p <- p + theme(panel.grid = element_blank())
p <- p + ylab(bquote(bold(R[PRS]^2)))
p <- p + theme(plot.title = element_text(size=7.5,family="sans",face="bold"), axis.text.x = element_text(family="sans",size=7.5,face="bold"), axis.text.y = element_text(family="sans",size=7.5), axis.title.x = element_text(family="sans",size=7.5,face = "bold"), axis.title.y = element_text(family="sans",size=7.5,face = "bold"))
p <- p + annotate("text",x = xpos, y= ypos,label=text,fontface=2,size=2,hjust=0)
p6 <- p

# Compare based on PRS P value
dat$P <- as.numeric(as.character(dat$P))
dat$PRS.P.adj <- p.adjust(dat$P, method = "fdr", n = length(dat$P))
dat$PRS.sig <- NA
dat[is.na(dat$PRS.P.adj) == FALSE & dat$PRS.P.adj < 0.05,"PRS.sig"] <- 1
dat[is.na(dat$PRS.P.adj) == FALSE & dat$PRS.P.adj >= 0.05,"PRS.sig"] <- 0
#dat[is.na(dat$PRS.P.adj) == FALSE & dat$P < 0.05,"PRS.sig"] <- 1
#dat[is.na(dat$PRS.P.adj) == FALSE & dat$P >= 0.05,"PRS.sig"] <- 0


# Fisher test

tryCatch({
testmat <- as.matrix(table(dat[,c("PRS.sig","JLIM")]))
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
test[test$PRS.sig == sigtit0 & test$JLIM == "-","Percentage"] <- round(100*(nrow(dat[dat$JLIMsimple == "JLIM-" & dat$PRS.sig == 0,])/nrow(dat[dat$JLIMsimple == "JLIM-",])),2)
test[test$PRS.sig == sigtit1 & test$JLIM == "-" ,"Percentage"] <- round(100*(nrow(dat[dat$JLIMsimple == "JLIM-" & dat$PRS.sig == 1,])/nrow(dat[dat$JLIMsimple == "JLIM-",])),2)
test[test$PRS.sig == sigtit0 & test$JLIM == "+","Percentage"] <- round(100*(nrow(dat[dat$JLIMsimple == "JLIM+" & dat$PRS.sig == 0,])/nrow(dat[dat$JLIMsimple == "JLIM+",])),2)
test[test$PRS.sig == sigtit1 & test$JLIM == "+" ,"Percentage"] <- round(100*(nrow(dat[dat$JLIMsimple == "JLIM+" & dat$PRS.sig == 1,])/nrow(dat[dat$JLIMsimple == "JLIM+",])),2)
jlimname0 <- as.character(dat[dat$JLIMsimple == "JLIM-","JLIM"][1])
jlimname1 <- as.character(dat[dat$JLIMsimple == "JLIM+","JLIM"][1])
test$JLIM <- c(jlimname0,jlimname1,jlimname0,jlimname1)
test$JLIM <- as.factor(test$JLIM)
test$PRS.sig <- as.factor(test$PRS.sig)
test$Percentage <- as.numeric(as.character((test$Percentage)))
#colnames(test)[2] <- PRS.pvalue

p <- ggplot(test, aes(x = JLIM, y = Percentage)) + geom_col(aes(fill = PRS.sig), width = 0.3,alpha=0.8)
p <- p + scale_color_jama() + scale_fill_jama(name="FDR PRS p",labels=c("<0.05",">=0.05"))
p <- p + xlab("") + ylab("%")
p <- p + theme(text = element_text(size=6,family="sans"), plot.title = element_text(size=5,family="sans",face="bold"), axis.text.x = element_text(family="sans",size=6,face="bold"), axis.text.y = element_text(family="sans",size=6), axis.title.x = element_text(family="sans",size=6,face = "bold"), axis.title.y = element_text(family="sans",size=6,face = "bold"))
p <- p + theme_minimal()
#p <- p + ggtitle(tit)
p <- p + theme(plot.title = element_text(size=7.5,family="sans",face="bold"), axis.text.x = element_text(family="sans",size=7.5,face="bold"), axis.text.y = element_text(family="sans",size=7.5), axis.title.x = element_text(family="sans",size=7.5,face = "bold"), axis.title.y = element_text(family="sans",size=7.5,face = "bold"))
p <- p + theme(legend.title=element_text(size=7.5,face="bold"),legend.text = element_text(size=7.5))
p <- p + theme(legend.key.size = unit(0.3, "cm"),legend.margin=margin(0,0,0,0),legend.box.margin=margin(0,-5,0,-30))
p <- p + theme(axis.title.y = element_text(vjust = -3))
p7 <- p

# Plot correlation JLIM minp to PRS P
dat$logPRSP <- -log10(dat$P)
#png(args[6],width=8,height=7,units="cm",res=300,type="cairo")
p <- ggplot(dat, aes(x=logminp,y=logPRSP))  + geom_point(color="grey50",alpha=0.65)
p <- p + xlab("-log(JLIM p)") + ylab("-log(PRS p)")
p <- p + geom_smooth(method=lm,color="black",level=0.95,size=0.5)
p <- p + theme(panel.background = element_rect(fill = "grey97"))
df <- dat[,c("logminp","logPRSP")]
colnames(df) <- c("x","y")
text <- as.character(paste(toString(lm_eqn(df,"log minp JLIM","log PRS p"))))
ypos <- max(dat$logPRSP,na.rm=T)*0.9
xpos <- max(dat$logminp,na.rm=T)*0.7
p <- p + theme_minimal()
p <- p + theme(plot.title = element_text(size=7.5,family="sans",face="bold"), axis.text.x = element_text(family="sans",size=6,face="bold"), axis.text.y = element_text(family="sans",size=6), axis.title.x = element_text(family="sans",size=6,face = "bold"), axis.title.y = element_text(family="sans",size=6,face = "bold"))
p <- p + annotate("text",x = xpos, y= ypos,label=text,fontface=2,size=2,hjust=0)
p8 <- p
#dev.off()





#png(args[9],width=150,height=120,units="mm",res=300,type="cairo")
png(args[9],width=180,height=120,units="mm",res=300,type="cairo")
grid.arrange(p1, p2, p3, p5, p6, p7, nrow = 2)
grid.text("a", x=unit(0.5, "mm"), y=unit(119, "mm"), just=c("left", "top"),gp=gpar(fontsize=7.5,fontface = "bold",family="sans"))
grid.text("b", x=unit(60.5, "mm"), y=unit(119, "mm"), just=c("left", "top"),gp=gpar(fontsize=7.5,fontface = "bold",family="sans"))
grid.text("c", x=unit(120.5, "mm"), y=unit(119, "mm"), just=c("left", "top"),gp=gpar(fontsize=7.5,fontface = "bold",family="sans"))
grid.text("d", x=unit(0.5, "mm"), y=unit(59, "mm"), just=c("left", "top"),gp=gpar(fontsize=7.5,fontface = "bold",family="sans"))
grid.text("e", x=unit(60.5, "mm"), y=unit(59, "mm"), just=c("left", "top"),gp=gpar(fontsize=7.5,fontface = "bold",family="sans"))
grid.text("f", x=unit(120.5, "mm"), y=unit(59, "mm"), just=c("left", "top"),gp=gpar(fontsize=7.5,fontface = "bold",family="sans"))
dev.off()









