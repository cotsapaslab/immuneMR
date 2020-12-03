args <- commandArgs(TRUE)

library(ggplot2)
library(gridExtra)
library(ggsci)
library(plyr)
library(grid)
library(plyr)
library(tidyr)
#--------------------------------
# TwoSampleMR
#--------------------------------

dat <- read.table(args[2],sep="\t",h=T)
jlim <- read.table(args[3],sep="\t",h=T)
method <- args[4]
pthres <- as.numeric(as.character(args[5]))
output <- args[6]

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
res <- merged[merged$method == method & as.numeric(as.character(merged$pthresMR)) == pthres,]
res$pval <- as.numeric(as.character(res$pval))
res$b <- as.numeric(as.character(res$b))
res$logP <- -log10(as.numeric(as.character(res$pval)))
res <- res[is.na(res$minp) == FALSE,]


# FDR pvalue
res$P.adj <- p.adjust(res$pval,method="fdr",n=length(res$pval))
res$sig <- 0
res[res$pval < 0.05,"sig"] <- 1
#res[res$P.adj < 0.05,"sig"] <- 1
print(table(res$sig))

# JLIM information
print(head(res))
pos <- res[res$JLIM == "JLIM+",]
neg <- res[res$JLIM == "JLIM-",]
print(nrow(pos))
print(nrow(neg))
res$JLIMname <- res$JLIM
res <- as.matrix(res)
res[res[,"JLIMname"] == "JLIM+","JLIMname"] <- paste0("shared effect\n(n=",as.character(nrow(pos)),")")
res[res[,"JLIMname"] == "JLIM-","JLIMname"] <- paste0("no shared effect\n(n=",nrow(neg),")")
res <- data.frame(res)
res$JLIMname <- as.factor(res$JLIMname)
res$pval <- as.numeric(as.character(res$pval))
res$b <- as.numeric(as.character(res$b))
res$logP <- as.numeric(as.character(res$logP))
res$absbeta <- abs(res$b)


# Fisher test and plot for significant TWMR results
dat <- res
test <- matrix(nrow=2,ncol=2)
test[1,] <- table(pos$sig)
test[2,] <- table(neg$sig)
print ("Fisher test number of sig TSMR results:")
print(test)
print(fisher.test(test))
pval <- fisher.test(test)$p.val
pval <- formatC(pval, format = "e", digits = 2)
tit <- paste0("Fisher's exact test p = ",pval)

test <- matrix(nrow=4,ncol=3,NA)
colnames(test) <- c("JLIM","TSMR.sig","Percentage")
test <- data.frame(test)
jlim0 <- paste0("no shared effect\n(n= ",nrow(dat[dat$JLIM == "JLIM-",]),")")
jlim1 <- paste0("shared effect\n(n= ",nrow(dat[dat$JLIM == "JLIM+",]),")")
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
p <- p + scale_color_jama() + scale_fill_jama(name="p",labels=c("<0.05",">=0.05"))
p <- p + xlab("") + ylab("%")
p <- p + theme_minimal()
p <- p + labs(fill = "p")
p <- p + theme(plot.title = element_text(size=7.5,family="sans",face="bold"), axis.text.x = element_text(family="sans",size=7.5,face="bold"), axis.text.y = element_text(family="sans",size=7.5), axis.title.x = element_text(family="sans",size=7.5,face = "bold"), axis.title.y = element_text(family="sans",size=7.5,face = "bold"))
p <- p + theme(legend.title=element_text(size=7.5,face="bold"),legend.text = element_text(size=7.5))
p <- p + theme(legend.key.size = unit(0.3, "cm"),legend.margin=margin(0,0,0,0),legend.box.margin=margin(0,-5,0,-30))
p <- p + theme(axis.title.y = element_text(vjust = -3))
p1 <- p











# Comparison of betas between JLIM+ and JLIM- combination
pval <- wilcox.test(res[res$JLIM == "JLIM+","b"],res[res$JLIM == "JLIM-","b"])$p.val
pval <- formatC(pval, format = "e", digits = 2)
tit <- paste0("M.-W.-Wilcoxon test p = ",as.character(pval))
print("beta:")
print(tit)

res$logabsbeta <- log10(res$absbeta)
p <- ggplot(res, aes(x=JLIMname, y=logabsbeta, color=JLIMname, fill=JLIMname)) + geom_boxplot(alpha=0.8,size=0.3)
p <- p + scale_color_nejm() + scale_fill_nejm()
p <- p + xlab("") + ylab(bquote(bold('Abs. est. causal effect (log'[10]~')')))
p <- p + ggtitle(tit)
p <- p + theme_minimal()
p <- p + theme(panel.grid = element_blank())
p <- p + theme(legend.title = element_blank())
p <- p + theme(legend.position = "none")
p <- p + theme(text = element_text(size=7.5,family="sans"), plot.title = element_text(size=7.5,family="sans",face="bold"), axis.text.x = element_text(family="sans",size=7.5,face="bold"), axis.text.y = element_text(family="sans",size=7.5), axis.title.x = element_text(family="sans",size=7.5,face = "bold"), axis.title.y = element_text(family="sans",size=7.5,face = "bold"))
p <- p + ggtitle(tit)
p2 <- p

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

p <- ggplot(res, aes(x=logminp,y=absbeta))  + geom_point(color="grey50",alpha=0.65,size=1)
p <- p + xlab(bquote(bold('-log'[10]~'(JLIM p)'))) + ylab("Abs. est. causal effect")
p <- p + geom_smooth(method=lm,color="black",level=0.95,size=0.5)
p <- p + theme(panel.background = element_rect(fill = "grey97"))
df <- res[,c("logminp","absbeta")]
colnames(df) <- c("x","y")
text <- as.character(paste(toString(lm_eqn(df,"log minp JLIM","beta"))))
ypos <- max(res$absbeta)*0.9
ypos2 <- max(res$absbeta)*0.82
xpos <- max(res$logminp)*0.7
p <- p + theme_minimal()
p <- p + theme(panel.grid = element_blank())
#
p <- p + annotate("text",x = xpos, y= ypos,label=text,size=2,fontface=2,hjust=0)
p <- p + theme(text = element_text(size=7.5,family="sans"), plot.title = element_text(size=7.5,family="sans",face="bold"), axis.text.x = element_text(family="sans",size=7.5,face="bold"), axis.text.y = element_text(family="sans",size=7.5), axis.title.x = element_text(family="sans",size=7.5,face = "bold"), axis.title.y = element_text(family="sans",size=7.5,face = "bold"))
p3 <- p


#--------------------------------
# TWMR
#--------------------------------

res <- read.table(args[1],sep="",h=T)
jlim <- read.table(args[3],sep="\t",h=T)

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

# JLIM information
dat$JLIMname <- dat$JLIM
dat <- as.matrix(dat)
dat[dat[,"JLIMname"] == "JLIM+","JLIMname"] <- paste0("shared effect\n(n=",as.character(nrow(pos)),")")
dat[dat[,"JLIMname"] == "JLIM-","JLIMname"] <- paste0("no shared effect\n(n=",nrow(neg),")")
dat <- data.frame(dat)
dat$JLIMname <- as.factor(dat$JLIMname)

dat$alpha <- as.numeric(as.character(dat$alpha))
dat$absalpha <- abs(dat$alpha)

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
jlim0 <- paste0("no shared effect","\n","(n= ",nrow(dat[dat$JLIM == "JLIM-",]),")")
jlim1 <- paste0("shared effect","\n","(n= ",nrow(dat[dat$JLIM == "JLIM+",]),")")
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
p <- p + scale_color_jama() + scale_fill_jama(name="FDR p",labels=c("<0.05",">=0.05"))
p <- p + xlab("") + ylab("%")
p <- p + theme_minimal()
p <- p + labs(fill = "FDR TWMR p")
#p <- p + ggtitle(tit)
p <- p + theme(plot.title = element_text(size=7.5,family="sans",face="bold"), axis.text.x = element_text(family="sans",size=7.5,face="bold"), axis.text.y = element_text(family="sans",size=7.5), axis.title.x = element_text(family="sans",size=7.5,face = "bold"), axis.title.y = element_text(family="sans",size=7.5,face = "bold"))
p <- p + theme(legend.title=element_text(size=7.5,face="bold"),legend.text = element_text(size=7.5))
p <- p + theme(legend.key.size = unit(0.3, "cm"),legend.margin=margin(0,0,0,0),legend.box.margin=margin(0,-5,0,-30))
p <- p + theme(axis.title.y = element_text(vjust = -3))
p4 <- p

# Comparison alpha
dat$alpha <- as.numeric(as.character(dat$alpha))
dat$absalpha <- abs(dat$alpha)

pval <- wilcox.test(dat[dat$JLIM == "JLIM+","absalpha"],dat[dat$JLIM == "JLIM-","absalpha"])$p.val
pval <- formatC(pval, format = "e", digits = 2)
tit <- paste0("M.-W.-Wilcoxon test p = ",as.character(pval))
print("Alpha:")
print(tit)


dat$logabsalpha <- log10(dat$absalpha)
p <- ggplot(dat, aes(x=JLIMname, y=logabsalpha, color=JLIMname, fill=JLIMname)) + geom_boxplot(trim=FALSE, alpha=0.8,size=0.3)
p <- p + scale_color_nejm() + scale_fill_nejm()
p <- p + xlab("") + ylab(bquote(bold('Abs. est. causal effect (log'[10]~')')))
p <- p + ggtitle(tit)
p <- p + theme_minimal()
p <- p + theme(legend.title = element_blank())
p <- p + theme(panel.grid = element_blank())
p <- p + theme(legend.position = "none")
p <- p + theme(text = element_text(size=7.5,family="sans"), plot.title = element_text(size=7.5,family="sans",face="bold"), axis.text.x = element_text(family="sans",size=7.5,face="bold"), axis.text.y = element_text(family="sans",size=7.5), axis.title.x = element_text(family="sans",size=7.5,face = "bold"), axis.title.y = element_text(family="sans",size=7.5,face = "bold"))

p5 <- p


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

p <- ggplot(dat, aes(x=logminp,y=absalpha))  + geom_point(color="grey50",alpha=0.65,size=1)
p <- p + xlab(bquote(bold('-log'[10]~'(JLIM p)'))) + ylab("Abs. est. causal effect")
p <- p + geom_smooth(method=lm,color="black",level=0.95,size=0.5)
p <- p + theme(panel.background = element_rect(fill = "grey97"))
df <- dat[,c("logminp","absalpha")]
colnames(df) <- c("x","y")
text <- as.character(paste(toString(lm_eqn(df,"log minp JLIM","absolute alpha"))))
p <- p + theme_minimal()
p <- p + theme(panel.grid = element_blank())
ypos <- max(dat$absalpha)*0.9
ypos2 <- max(dat$absalpha)*0.82
xpos <- max(dat$logminp)*0.7
p <- p + annotate("text",x = xpos, y= ypos,label=text,size=2,fontface=2,hjust=0)
p <- p + theme(text = element_text(size=7.5,family="sans"), plot.title = element_text(size=7.5,family="sans",face="bold"), axis.text.x = element_text(family="sans",size=7.5,face="bold"), axis.text.y = element_text(family="sans",size=7.5), axis.title.x = element_text(family="sans",size=7.5,face = "bold"), axis.title.y = element_text(family="sans",size=7.5,face = "bold"))
p6 <- p



#png(output,height=100,width=120,units="mm",res=300,type="cairo")
png(output,height=100,width=180,units="mm",res=300,type="cairo")
grid.arrange(p2, p3, p1, p5, p6, p4, nrow = 2)
#grid.arrange(p2, p3, p5, p6, nrow = 2)
grid.text("a", x=unit(0, "mm"), y=unit(99, "mm"), just=c("left", "top"),gp=gpar(fontsize=7.5,fontface = "bold",family="sans"))
grid.text("b", x=unit(60, "mm"), y=unit(99, "mm"), just=c("left", "top"),gp=gpar(fontsize=7.5,fontface = "bold",family="sans"))
grid.text("c", x=unit(120, "mm"), y=unit(99, "mm"), just=c("left", "top"),gp=gpar(fontsize=7.5,fontface = "bold",family="sans"))
grid.text("d", x=unit(0, "mm"), y=unit(49, "mm"), just=c("left", "top"),gp=gpar(fontsize=7.5,fontface = "bold",family="sans"))
grid.text("e", x=unit(60, "mm"), y=unit(49, "mm"), just=c("left", "top"),gp=gpar(fontsize=7.5,fontface = "bold",family="sans"))
grid.text("f", x=unit(120, "mm"), y=unit(49, "mm"), just=c("left", "top"),gp=gpar(fontsize=7.5,fontface = "bold",family="sans"))
dev.off()
