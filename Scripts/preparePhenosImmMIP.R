args <- commandArgs(TRUE)

library(readxl)
library(ggplot2)
library(grid)
library(knitr)
library(png)
library(plyr)
library(methods)

datafolder <- args[1]
wd <- args[1]
source(paste0("/home/cg859/scratch60/Scripts/normalityAnalysis.r",collapse=""))

#functions 
"rntransform" <-
		function(formula,data,family=gaussian) {
	
	if ( is(try(formula,silent=TRUE),"try-error") ) { 
		if ( is(data,"gwaa.data") ) data1 <- phdata(data)
		else if ( is(data,"data.frame") ) data1 <- data
		else stop("'data' must have 'gwaa.data' or 'data.frame' class")
		formula <- data1[[as(match.call()[["formula"]],"character")]] 
	}
	
	var <- ztransform(formula,data,family)
	out <- rank(var) - 0.5
	out[is.na(var)] <- NA
	mP <- .5/max(out,na.rm=T)
	out <- out/(max(out,na.rm=T)+.5)
	out <- qnorm(out)
	out
}

"ztransform" <- 
function(formula,data,family=gaussian) {
	if (missing(data)) {
		if(is(formula,"formula")) 
			data <- environment(formula)
		else  
			data <- environment()
#		wasdata <- 0
	} else {
		if (is(data,"gwaa.data")) {
			data <- data@phdata
		} 
		else if (!is(data,"data.frame")) {
			stop("data argument should be of gwaa.data or data.frame class")
		}
#		attach(data,pos=2,warn.conflicts=FALSE)
#		wasdata <- 1
	}
	
	if (is.character(family)) 
           family <- get(family, mode = "function", envir = parent.frame())
	if (is.function(family)) 
           family <- family()
	if (is.null(family$family)) {
           print(family)
           stop("'family' not recognized")
	}
	
	if ( is(try(formula,silent=TRUE),"try-error") ) { 
		formula <- data[[as(match.call()[["formula"]],"character")]] 
	}
	
	if (is(formula,"formula")) {
#		mf <- model.frame(formula,data,na.action=na.omit,drop.unused.levels=TRUE)
		mf <- model.frame(formula,data,na.action=na.pass,drop.unused.levels=TRUE)
		mids <- complete.cases(mf)
		mf <- mf[mids,]
		y <- model.response(mf)
		desmat <- model.matrix(formula,mf)
		lmf <- glm.fit(desmat,y,family=family)
#		if (wasdata) 
#			mids <- rownames(data) %in% rownames(mf)
#		else 
		resid <- lmf$resid
#		print(formula)
	} else if (is(formula,"numeric") || is(formula,"integer") || is(formula,"double")) {
		y <- formula
		mids <- (!is.na(y))
		y <- y[mids]
		resid <- y
		if (length(unique(resid))==1) stop("trait is monomorphic")
		if (length(unique(resid))==2) stop("trait is binary")
	} else {
		stop("formula argument must be a formula or one of (numeric, integer, double)")
	}
	y <- (resid-mean(resid))/sd(resid)
#	if (wasdata==1) detach(data)
	tmeas <- as.logical(mids)
	out <- rep(NA,length(mids))
	out[tmeas] <- y
	out
}





phen <- read.table(paste0(datafolder,"/Pheno/ImmPhen/LabExMI_rawfacs.txt",collapse=""),sep="\t",header=T)
raw <- phen 

first <- 2
last <- ncol(phen)
raw <- phen 


# Rank transformation

for (i in first:last)
  {
  phen[paste0(colnames(phen)[i],"_norm")] <- NA
  phen[paste0(colnames(phen)[i],"_norm")] <- rntransform(phen[,i])
  }

pheno <- phen
pheno <- pheno[,168:ncol(pheno)]
first <- 1
last <- ncol(pheno)

# Check normality of transformed distributions
notnorm <- cbind(colnames(pheno),colnames(pheno))
notnorm[,2] <- NA
for (i in 1:nrow(notnorm)){
val <- notnorm[i,1]
notnorm[i,2] <- shapiro.test(as.numeric(as.character(pheno[,val])))$p.value}
notnorm <- notnorm[as.numeric(as.character(notnorm[,2])) < 0.0003012048,]


# Remove two traits based on plots
pheno$N_HLADRposCD56hi.panel4_norm <- NULL
pheno$N_CD4negCD8neg_NKT_HLADRpos.panel3_norm <- NULL

# Create binary phenotypes for two parameters (see Scripts/PRS/binarize_phenotypes.R for desicion plots)
traits <- c("N_FounderB.panel6","N_HLADRpos_in_CD4pos_EMRA.panel1")
for (i in 1:length(traits)){
trait <- traits[i]
test <- raw[,trait]
testq <- qqnorm(test)
x <- testq$x
y <- testq$y
out <- data.frame(x=x,y=y)
out <- out[is.na(out$y) == FALSE,]
out$s <- NA
    for (i in 1:(nrow(out)-1)){
    dx <- out[i+1,1]-out[i,1]
    dy <- out[i+1,2]-out[i,2]
    a <- dy/dx
    out[i,3] <- a
    }
out$ss <- NA
  for (i in 1:(nrow(out)-1)){
  dx <- out[i+1,1]-out[i,1]
  ds <- out[i+1,3]-out[i,3]
  a <- dy/ds
  out[i,4] <- a
}
out <- out[is.infinite(out$ss) == FALSE,]
pox = as.numeric(as.character(out[out$ss == max(out$ss,na.rm=T) & is.na(out$ss) == FALSE,"x"]))
out$x <- as.numeric(as.character(out$x))
poy = as.numeric(as.character(out[out$x == pox,"y"]))
newname <- paste0(trait,"_binary")
rntname <- paste0(trait,"_norm")
pheno <- data.frame(pheno)
pheno[,rntname] <- raw[,trait]
pheno[,newname] <- NA
pheno[is.na(pheno[,rntname]) == FALSE & pheno[,rntname] < poy,newname] <- 1
pheno[is.na(pheno[,rntname]) == FALSE & pheno[,rntname] >= poy,newname] <- 2
pheno[,rntname] <- NULL
}


pheno_trans <- cbind(raw[,1],pheno)
colnames(pheno_trans)[1] <- "IID"
pheno_trans[,"FID"] <- "LabExMI"
pheno_trans <- pheno_trans[,c(ncol(pheno_trans),1:(ncol(pheno_trans)-1))]


write.table(pheno_trans,paste0(datafolder,"/Pheno/ImmPhen/phenotypes_transformed.txt",collapse=""),sep="\t",row.names=F,quote=F)
write.table(colnames(pheno_trans[3:ncol(pheno_trans)]),paste0(datafolder,"/Pheno/ImmPhen/PhenoList.txt"),sep="\n",row.names=F,quote=F,col.names=F)

