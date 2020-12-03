# ===============================================
# Generation of plots for assessment of normality
# version 2019-03-04
# ===============================================
# Contains the following functions:
# ---------------------------------
# normality
# tablePlotMatrix (called by normality)

require(ggplot2)
require(grid)
require(reshape2)
require(RColorBrewer)

normality <- function(phenoData, nRows, nCols, first=1, last=ncol(phenoData), save=T, prefix="", mLog=F, mSqrt=F, pow=NA, titles="")
  {
  
  if(toString(titles)!="") colnames(phenoData)[first:last] <- titles

  if (mLog) # convert to log()
    {
    for (i in first:last){
      if (min(phenoData[,i],na.rm=T) > 0){phenoData[,i] <- log(phenoData[,i])}
      else{phenoData[,i] <- NA}   

    }
    phenoData <- phenoData[,colSums(is.na(phenoData))<nrow(phenoData)]
    last <- ncol(phenoData)
    outFile <- paste0(prefix,"normality_log.pdf")
    }
  else if (mSqrt) # convert to sqrt()
    {
    #phenoData[,first:last] <- sqrt(phenoData[,first:last])
    for (i in first:last){
      if (min(phenoData[,i],na.rm=T) >= 0){phenoData[,i] <- sqrt(phenoData[,i])}
      else{phenoData[,i] <- NA}
    }
    phenoData <- phenoData[,colSums(is.na(phenoData))<nrow(phenoData)]
    last <- ncol(phenoData)
    outFile <- paste0(prefix,"normality_sqrt.pdf")
    }
  else if (!is.na(pow)) # convert to x^pow
    {
    phenoData[,first:last] <- phenoData[,first:last]^pow
    phenoData <- phenoData[,colSums(is.na(phenoData))<nrow(phenoData)]
    last <- ncol(phenoData)
    outFile <- paste0(prefix,"normality_power",pow,".pdf")
    }
  else outFile <- paste0(prefix,"normality_standard.png") # plot original values

  if(save) png(outFile, width=2500, height=3500, units="px",pointsize=12,res=300,type="cairo")



  shap1 <- data.frame(sapply(phenoData[,first:last], function(x){t(unlist(shapiro.test(x)[1:2]))})) # apply shapiro test
  rownames(shap1) <- c("W","p")

  gp1 <- tablePlotMatrix(shap1, nCols=nCols, text=F, textSize=8)

  phenoData <- melt(phenoData,measure.vars=colnames(phenoData)[first:last])
  gp2 <- ggplot(phenoData,aes(x=value)) +
    geom_histogram(aes(y=..density..), fill=brewer.pal(12,"Paired")[1], colour=brewer.pal(12,"Paired")[2]) +
    geom_density(colour="grey40") +scale_y_continuous(name="Frequency") +
    scale_x_continuous(name="Measurement") +theme_bw() +
    theme(panel.grid.major=element_line(colour="grey60"),panel.grid.minor=element_blank()) +
    facet_wrap(~variable,scale="free",nrow=nRows,ncol=nCols)

  pushViewport(viewport(layout=grid.layout(nRows,nCols)))
  vplayout <- function(x,y) viewport(layout.pos.row=x, layout.pos.col=y)
  print(gp2, vp=vplayout(1:nRows, 1:nCols))
  print(gp1, vp=vplayout(nRows, nCols))
  if (save) dev.off()

  return(shap1)
  }

tablePlotMatrix <- function(data, nCols, text=T, textSize=4, bc=F) # data = shapiro.test
  {
  p <- as.data.frame(data[2,])
  if (bc) {bonf = 0.05/length(p)} else {bonf=0.05}

  #convert to symbols
  ps <- p
  ps[which(p>=bonf)] <- "ns"
  ps[which(p<bonf)] <- "*"
  ps[which(p<bonf/5)] <- "**"
  ps[which(p<bonf/50)] <- "***"

  ps2 <- as.data.frame(t(ps))
  ps2 <- cbind(ps2,as.data.frame(t(p)))
  ps2 <- cbind(ps2,rownames(ps2))

  n <- ncol(data)
  nRows <- ceiling(n/nCols)
  fullRows <- floor(n/nCols)
  lastRow <- n-(nCols*floor(n/nCols))
  c1 <- NULL

  for(i in 0:(fullRows-1)) {c1 <- c(c1, rep(nRows-i,nCols))}

  c1 <- c(c1, rep(1,lastRow))
  if (lastRow==0)
    { c2 <- c(rep(c(1:nCols),fullRows))
    } else c2 <- c(rep(c(1:nCols),fullRows), c(1:lastRow))
  ps2 <- cbind(ps2,c1,c2)
  colnames(ps2) <- c("p","num","name","row", "col")
  ps2$p <- factor(ps2$p,levels=c("*","**","***","ns"))

  pal1 <- c(brewer.pal(9,"Purples")[6:8],brewer.pal(9,"Set1")[3])
  if (text)
    { mPlot <- ggplot(ps2, aes(x=col, y=row)) +geom_tile(aes(fill=p),color="white") +geom_text(aes(label=name),color="white") +scale_fill_manual(values=pal1,breaks=c("*","**","***","ns"),drop=F,guide=F) +scale_x_continuous(c(1:4),expand=c(0,0)) +scale_y_continuous(c(1:4),expand=c(0,0)) +theme_bw() +theme(legend.title=element_blank(),axis.ticks=element_blank(),axis.text=element_blank(),panel.grid.minor=element_blank(),panel.grid.major=element_blank(),axis.title=element_blank())
    } else mPlot <- ggplot(ps2, aes(x=col, y=row)) +geom_tile(aes(fill=p),color="white") +geom_text(aes(label=p),color="white",size=textSize) +scale_fill_manual(values=pal1,breaks=c("*","**","***","ns"),drop=F,guide=F) +scale_x_continuous(c(1:4),expand=c(0,0)) +scale_y_continuous(c(1:4),expand=c(0,0))  +theme_bw() +theme(legend.title=element_blank(),axis.ticks=element_blank(),axis.text=element_blank(),panel.grid.minor=element_blank(),panel.grid.major=element_blank(),axis.title=element_blank())

  return(mPlot)
  }
