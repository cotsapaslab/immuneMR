args <- commandArgs(TRUE)

library(ggplot2)
library(genetics)
library(RColorBrewer)
library(gridExtra)
library(gCMAP)
library(trio)
library(grid)
library(lemon)

ass1file <- args[1]
ass2file <- args[2]
pheno <- args[3]
ldfile <- args[4]
ldplinkfile <- args[5]
leadfile <- args[6]
firsttrait <- args[7]
outfolder <- args[8]

pedfileG <- args[9]
mapfileG <- args[10]
phenofileG <- args[11]
pedfileI <- args[12]
mapfileI <- args[13]
phenofileI <- args[14]
firsttraitn <- args[15]
firsttraitp <- args[16]
DS <- args[17]

ass1 <- read.table(ass1file,sep="",header=T)
ass2 <- read.table(toString(ass2file),sep="",header=T)
ld <- read.table(ldfile,sep="")
phase <- read.table(ldplinkfile,sep="",header=T)
leadf <- read.table(leadfile,sep="")
lead1 <- toString(leadf[1,1])
chrom <- toString(ass1[1,"CHR"])
print("hallo1")
pedG <- read.pedfile(pedfileG)
mapG <- read.table(mapfileG)
phenoG <- read.table(phenofileG,sep="\t",h=T)
print("hallo2")
print(pedfileI)
pedI <- read.pedfile(pedfileI)
print("hallo2.2")
mapI <- read.table(mapfileI)
phenoI <- read.table(phenofileI,sep="\t",h=T)
print("hallo3")

if ("BETA" %in% colnames(ass1) == TRUE && "SE" %in% colnames(ass1) == TRUE){
ass1$Z.x <- ass1$BETA/ass1$SE
}else if ( "OR" %in% colnames(ass1) == TRUE && "SE" %in% colnames(ass1) == TRUE){
ass1$Z.x <- log(ass1$OR)/ass1$SE 
}else{
ass1$Z.x <- zScores(ass1$P) 
}

ass1[is.na(ass1$OR) == FALSE & ass1$OR < 1,"Z.x"] <- ass1[is.na(ass1$OR) == FALSE & ass1$OR < 1,"Z.x"] * -1


if ("TEST" %in% colnames(ass2) == TRUE){
ass2 <- ass2[ass2[,"TEST"] == "ADD",]}

dat <- merge(ass1,ass2,by="SNP")
dat$logPx <- -log10(dat$P.x)
dat$logPy <- -log10(dat$P.y)

if ("BETA.y" %in% colnames(dat) == TRUE && "SE.y" %in% colnames(dat)){
dat$Z.y <- dat$BETA.y/dat$SE.y
}else if ("OR.y" %in% colnames(dat) == TRUE){
dat$Z.y <- log(dat$OR.y)/dat$SE.y
}else{
dat$Z.y <- zScores(dat$P.y)
}



dat$BP_MB <- dat$BP.x/1000000


dat <- dat[dat$BP.x %in% ld[,2] == TRUE,]
#dat <- dat[dat[,1] %in% ld[,3] == TRUE | dat[,3] %in% ld[,1] == TRUE,]

#calculate ld to MS lead SNP
ld <- ld[ld[,2] %in% dat$BP.x == T,]
dat$r2tolead <- NA
lead <- lead1
leadpos <- dat[dat$SNP == lead,"BP.x"]
if (lead %in% dat$SNP == FALSE){
lead <- dat[dat$P.x == min(dat$P.x),"SNP"]	
leadpos <- dat[dat$SNP == lead,"BP.x"]
}	
	snp <- toString(lead)
	i <- 8
	gen <- c()
	while (i <= ncol(ld)){
	#genn <- paste0(toString(ld[ld[,3] == toString(snp),i]),"/",toString(ld[ld[,3] == toString(snp),i+1]))
	genn <- paste0(toString(ld[ld[,2] == toString(leadpos),i]),"/",toString(ld[ld[,2] == toString(leadpos),i+1]))
	gen <- c(gen,genn)
	i <- i + 2
	}
	glead <- genotype(gen)


for (j in 1:nrow(dat)){
snp <- toString(dat[j,"SNP"])
snppos <- toString(dat[j,"BP.x"])
ld[grepl(snp,ld[,3]) == TRUE,3] <- snp
	i <- 8
        gen <- c()
	while (i <= ncol(ld)){
        #genn <- paste0(toString(ld[ld[,3] == toString(snp),i]),"/",toString(ld[ld[,3] == toString(snp),i+1]))
        genn <- paste0(toString(ld[ld[,2] == toString(snppos),i]),"/",toString(ld[ld[,2] == toString(snppos),i+1]))
	gen <- c(gen,genn)
        i <- i + 2
 	}
	g1 <- NA
	g1 <- genotype(gen)
	r2 <- NA
	tryCatch({
	r2 <- LD(g1,glead)$r * LD(g1,glead)$r
	},error=function(cond){})
	dat[j,"r2tolead"] <- r2	

}

# flip to MS alleles
dat$A1.y <- as.character(dat$A1.y)
dat$A1.x <- as.character(dat$A1.x)
dat$A2 <- as.character(dat$A2)

dat[dat[,"A1.y"] != dat[,"A1.x"] & dat[,"A1.y"] == dat[,"A2"],"Z.y"] <- dat[dat[,"A1.y"] != dat[,"A1.x"] & dat[,"A1.y"] == dat[,"A2"],"Z.y"] * -1
dat[dat[,"A1.y"] != dat[,"A1.x"] & dat[,"A1.y"] != dat[,"A2"],"Z.y"] <- NA
dat <- dat[is.na(dat[,"Z.x"]) == FALSE & is.na(dat[,"Z.y"]) == FALSE,]

#flip to inphase-alleles
 phase$A <- NA
 phase$B <- NA
 phase$C <- NA
 phase$D <- NA
 dat$flip <- NA
 phase <- as.matrix(phase)
 phase <- phase[phase[,"SNP_A"] == lead1 | phase[,"SNP_B"] == lead1 ,]
 phase[,"A"] <- substr(phase[,"PHASE"],1,1)
 phase[,"B"] <- substr(phase[,"PHASE"],2,2)
 phase[,"C"] <- substr(phase[,"PHASE"],4,4)
 phase[,"D"] <- substr(phase[,"PHASE"],5,5)
 test <- phase[phase[,"SNP_A"] == lead1,]
 pall <- test[1,"A"]
 if (dat[dat[,"SNP"] == lead1,"A1.x"] == pall){
	checkall <- "eins"
 }else{
        checkall <- "zwei"
}
check <- checkall

phase <- data.frame(phase)
for (i in 1:nrow(dat)){
	snp <- dat[i,"SNP"]
	if (snp != lead1){
		test <- phase[phase[,"SNP_A"] == toString(lead1) & phase[,"SNP_B"] == toString(snp),]
	if (check == "eins" && test[,"B"] == dat[dat[,"SNP"] == snp,"A1.x"] && test[,"D"] == dat[dat[,"SNP"] == snp,"A2"]){dat[i,"flip"] <- "no"}
	if (check == "eins" && test[,"B"] == dat[dat[,"SNP"] == snp,"A2"] && test[,"D"] == dat[dat[,"SNP"] == snp,"A1.x"]){dat[i,"flip"] <- "yes"}
	if (check == "zwei" && test[,"B"] == dat[dat[,"SNP"] == snp,"A1.x"] && test[,"D"] == dat[dat[,"SNP"] == snp,"A2"]){dat[i,"flip"] <- "yes"}
        if (check == "zwei" && test[,"B"] == dat[dat[,"SNP"] == snp,"A2"] && test[,"D"] == dat[dat[,"SNP"] == snp,"A1.x"]){dat[i,"flip"] <- "no"}
}}	

dat[dat[,"flip"] == "yes" & is.na(dat[,"flip"]) == FALSE & is.na(dat[,"Z.x"]) == FALSE,"Z.x"] <- dat[dat[,"flip"] == "yes" & is.na(dat[,"flip"]) == FALSE & is.na(dat[,"Z.x"]) == FALSE,"Z.x"]*-1
dat[dat[,"flip"] == "yes" & is.na(dat[,"flip"]) == FALSE & is.na(dat[,"Z.y"]) == FALSE,"Z.y"] <- dat[dat[,"flip"] == "yes" & is.na(dat[,"flip"]) == FALSE & is.na(dat[,"Z.y"]) == FALSE,"Z.y"]*-1
dat[is.na(dat[,"flip"]) == TRUE & is.na(dat[,"Z.y"]) == FALSE,"Z.y"] <- NA
dat[is.na(dat[,"flip"]) == TRUE & is.na(dat[,"Z.x"]) == FALSE,"Z.x"] <- NA



Gene <- pheno
if (grepl("ENSG",pheno) == TRUE){
        library("biomaRt")
        grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
        ens_name <- unlist(strsplit(pheno,"\\."))[1]
        ann <- getBM(c("hgnc_symbol","description","chromosome_name","band","strand","start_position","end_position","ensembl_gene_id"), "ensembl_gene_id", ens_name, grch37)
        if (nrow(ann) > 0){pheno <- ann[1,"hgnc_symbol"]}

}

pheno <- gsub("._norm","",pheno)
pheno <- gsub("_norm","",pheno)
firsttrait <- unlist(strsplit(firsttrait,".panel"))[1]

# Prepare data for genotype-phenotype boxplots
ped <- pedG; map <- mapG
rownames(ped) <- as.character(ped[,2])
ped <- ped[,7:(ncol(ped))]
j <- ncol(ped)+1
i <- 1
max <- ncol(ped)
while(i < max){
for (n in 1:nrow(ped)){
ped[n,j] <- paste0(as.character(ped[n,i]),as.character(ped[n,i+1]))}
i <- i + 2
j <- j + 1}
ped <- ped[,(max+1):ncol(ped)]
colnames(ped) <- as.character(map[,2])
ped$lead <- ped[,lead]
ped <- ped[ped$lead != "00" & is.na(ped$lead) == FALSE,]
ped$IID <- rownames(ped)
ped <- ped[,colnames(ped) == "IID" | colnames(ped) == lead,]
ped <- merge(ped,phenoG[,c("IID",Gene)],by="IID")
colnames(ped)  <- c("IID","lead","pheno")
genophenoG <- ped
print(head(genophenoG))


ped <- pedI; map <- mapI
rownames(ped) <- as.character(ped[,2])
ped <- ped[,7:(ncol(ped))]
j <- ncol(ped)+1
i <- 1
max <- ncol(ped)
while(i < max){
for (n in 1:nrow(ped)){
ped[n,j] <- paste0(as.character(ped[n,i]),as.character(ped[n,i+1]))}
i <- i + 2
j <- j + 1}
ped <- ped[,(max+1):ncol(ped)]
colnames(ped) <- as.character(map[,2])
ped$lead <- ped[,lead]
ped <- ped[ped$lead != "00" & is.na(ped$lead) == FALSE,]
ped$IID <- rownames(ped)
ped <- ped[,colnames(ped) == "IID" | colnames(ped) == lead,]
ped <- merge(ped,phenoI[,c("IID",firsttraitn)],by="IID")
colnames(ped)  <- c("IID","lead","pheno")
genophenoI <- ped



# Man1 Plot 
pman1 <- ggplot(dat,aes(x=BP_MB,y=logPx,color=r2tolead)) + geom_point(size=2,alpha=0.65,shape=16) 
pman1 <- pman1 + scale_color_gradient(low="grey50", high="red",name="")
pman1 <- pman1 + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
pman1 <- pman1 + ylab(bquote(bold('-log'[10]~'(p)'))) + xlab("")
pman1 <- pman1 + theme_minimal()
pman1 <- pman1 + theme(axis.line=element_line(), axis.ticks = element_line())
pman1 <- pman1 + coord_capped_cart(bottom='both', left='both')
pman1 <- pman1 + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
pman1 <- pman1 + theme(axis.line.x = element_line(color="white"))
pman1 <- pman1 + theme(panel.grid = element_blank())
pman1 <- pman1 + theme(axis.title.y =  element_text(family="sans",size=8,face="bold")) 
png(paste0(outfolder,"Figure_manhattan_",firsttrait,"_withlegend.png"),width=100,height=60,units="mm",type="cairo",res=300)
pman1
dev.off()
pman1 <- pman1 + theme(legend.position = "none") 
png(paste0(outfolder,"Figure_manhattan_",firsttrait,"_nolegend.png"),width=100,height=60,units="mm",type="cairo",res=300)
pman1
dev.off()

# Man2 Plot
pman2 <- ggplot(dat,aes(x=BP_MB,y=logPy,color=r2tolead)) + geom_point(size=2,alpha=0.65,shape=16) 
pman2 <- pman2 + scale_color_gradient(low="grey50", high="red",name="")
pman2 <- pman2 + ylab(bquote(bold('-log'[10]~'(eQTL p)')))
pman2 <- pman2 + theme_minimal()
pman2 <- pman2 + theme(axis.line=element_line(), axis.ticks = element_line())
pman2 <- pman2 + coord_capped_cart(bottom='both', left='both')
pman2 <- pman2 + theme(panel.grid = element_blank())
pman2 <- pman2 + theme(axis.title.y =  element_text(family="sans",size=8,face="bold"), axis.text.y = element_text(family="sans",size=8)) 
pman2 <- pman2 + theme(legend.position = "none") 

pman21 <- pman2 + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
pman21 <- pman21 + theme(axis.line.x = element_line(color="white"))

#pman2 <- pman2 + scale_x_continuous(breaks=c(169.60,169.75))
#pman2 <- pman2 + scale_x_continuous(breaks=c(207.55,207.75))
pman2 <- pman2 + xlab(paste0("Chromosome ",chrom," position"))
pman2 <- pman2 + theme(axis.title.x =  element_text(family="sans",size=8,face="bold"), axis.text.x = element_text(family="sans",size=8) ) 

png(paste0(outfolder,"Figure_manhattan_",firsttrait,"_",pheno,"_",DS,".png"),width=100,height=60,units="mm",type="cairo",res=300)
pman21
dev.off()

png(paste0(outfolder,"Figure_manhattan_",firsttrait,"_",pheno,"_",DS,"_Xaxis.png"),width=100,height=60,units="mm",type="cairo",res=300)
pman2
dev.off()

library(grid)
gpman1 <- ggplot_build(pman1)
gpman2 <- ggplot_build(pman2)
gpman21 <- ggplot_build(pman21)
ggpman1 <- ggplot_gtable(gpman1)
ggpman2 <- ggplot_gtable(gpman2)
ggpman21 <- ggplot_gtable(gpman21)
g1 <- gtable:::rbind_gtable(ggpman1,ggpman2,"last") #last bei CR2 plot
g2 <- gtable:::rbind_gtable(ggpman1,ggpman21,"last") #last bei CR2 plot

outname1 <- paste0(outfolder,"Figure_manhattan_",pheno,"_",DS,".png")
outname2 <- paste0(outfolder,"Figure_manhattan_",pheno,"_",DS,"_Xaxis.png")

png(paste0(outname1),width=100,height=120,units="mm",type="cairo",res=300)
grid.newpage()
grid.draw(g1)
dev.off()
png(paste0(outname2),width=100,height=120,units="mm",type="cairo",res=300)
grid.newpage()
grid.draw(g2)
dev.off()



# Z Plot 
pz <- ggplot(dat,aes(x=Z.y,y=Z.x)) + geom_point(col="black",size=2,alpha=0.55,shape=16)
pz <- pz + xlab("Z Gene exp.")
pz <- pz + theme(axis.title.y =  element_text(family="sans",size=8,face="bold")) + scale_color_gradient(low="grey50", high="red",name="")
pz <- pz + theme_minimal()
pz <- pz + theme(axis.line=element_line(), axis.ticks = element_line())
pz <- pz + coord_capped_cart(bottom='both', left='both')
pz <- pz + theme(legend.position = "none")
pz <- pz + theme(panel.grid = element_blank())
pz <- pz + theme(axis.title.y =  element_text(family="sans",size=8,face="bold"), axis.text.y = element_text(family="sans",size=8))
pz <- pz + ylab("Z ImmPhen")
pz1 <- pz + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
pz1 <- pz1 + theme(axis.line.x = element_line(color="white"))

#pz <- pz + scale_x_continuous(breaks=c(-4,-2,0))
#pz <- pz + xlab("Z (L-selectin in neutrophils)")
pz <- pz + theme(axis.title.x =  element_text(family="sans",size=8,face="bold"), axis.text.x = element_text(family="sans",size=8) ) 

gpz1 <- ggplot_build(pz1)
gpz <- ggplot_build(pz)
ggpz1 <- ggplot_gtable(gpz1)
ggpz <- ggplot_gtable(gpz)
gz1 <- gtable:::rbind_gtable(ggpman2,ggpz1,"first")
gz <- gtable:::rbind_gtable(ggpman2,ggpz,"first")
outname3 <- paste0(outfolder,"Figure_Z_",firsttrait,"_",pheno,"_",DS,".png")
outname4 <- paste0(outfolder,"Figure_Z_",firsttrait,"_",pheno,"_",DS,"_Xaxis.png")
png(paste0(outname3),width=50,height=120,units="mm",type="cairo",res=300)
grid.newpage()
grid.draw(gz)
dev.off()
png(paste0(outname4),width=50,height=120,units="mm",type="cairo",res=300)
grid.newpage()
grid.draw(gz1)
dev.off()


png(paste0(outfolder,"Figure_Z_",firsttrait,"_",pheno,"_",DS,".png"),width=50,height=60,units="mm",type="cairo",res=300)
pz1
dev.off()

png(paste0(outfolder,"Figure_Z_",firsttrait,"_",pheno,"_",DS,"_Xaxis.png"),width=50,height=60,units="mm",type="cairo",res=300)
pz
dev.off()


# Box Plot 1
bp1 <- ggplot(genophenoI,aes(x=lead,y=pheno,fill=lead)) + geom_boxplot() + scale_fill_grey() 
bp1 <- bp1 + ylab("rank transformed ImmPhen")
bp1 <- bp1 + theme_minimal()
bp1 <- bp1 + theme(legend.position = "none")
bp1 <- bp1 + theme(axis.line=element_line(), axis.ticks = element_line())
bp1 <- bp1 + coord_capped_cart(bottom='both', left='both')
bp1 <- bp1 + theme(panel.grid = element_blank())
bp1 <- bp1 + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
bp1 <- bp1 + theme(axis.line.x = element_line(color="white"))
bp1 <- bp1 + theme(axis.title.y =  element_text(family="sans",size=8,face="bold"), axis.text.y =  element_text(family="sans",size=8,face="bold"))
png(paste0(outfolder,"Figure_Boxplot_",firsttrait,".png"),width=50,height=60,units="mm",type="cairo",res=300)
bp1
dev.off()


# Box Plot 2
bp2 <- ggplot(genophenoG,aes(x=lead,y=pheno,fill=lead)) + geom_boxplot() + scale_fill_grey() 
bp2 <- bp2 + ylab("rank transformed expression values")
bp2 <- bp2 + theme_minimal()
bp2 <- bp2 + theme(legend.position = "none")
bp2 <- bp2 + theme(axis.line=element_line(), axis.ticks = element_line())
bp2 <- bp2 + coord_capped_cart(bottom='both', left='both')
bp2 <- bp2 + theme(panel.grid = element_blank())
#bp2 <- bp2 + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
bp2 <- bp2 + theme(axis.title.y =  element_text(family="sans",size=8,face="bold"), axis.text.y =  element_text(family="sans",size=8,face="bold")) 

bp21 <- bp2 + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
bp21 <- bp21 + theme(axis.line.x = element_line(color="white"))

bp2 <- bp2 + xlab(lead)
bp2 <- bp2 + theme(axis.title.x =  element_text(family="sans",size=8,face="bold"), axis.text.x =  element_text(family="sans",size=8,face="bold")) 

png(paste0(outfolder,"Figure_Boxplot_",firsttrait,"_",pheno,"_",DS,".png"),width=50,height=60,units="mm",type="cairo",res=300)
bp21
dev.off()

png(paste0(outfolder,"Figure_Boxplot_",firsttrait,"_",pheno,"_",DS,"_Xaxis.png"),width=50,height=60,units="mm",type="cairo",res=300)
bp2
dev.off()

