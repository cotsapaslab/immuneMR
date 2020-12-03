args <- commandArgs(TRUE)

#module load R/3.5.0-foss-2016b-avx2

exposurefilename <- args[1]
outcomefilename <- args[2]
genexpds <- args[3]
jlim <- args[4]
outputfolder <- args[5]

library(gridExtra)
library(TwoSampleMR)
library(ggplot2)

# Read exposure data file
dat <- read.table(exposurefilename,sep="",h=T)
keep <- dat

pthreslist <- c(0.00000005,0.0000001,0.0000005,0.000001,0.000005,0.00001,0.00005,0.0001)
firstres <- 1
firstplei <- 1
firsthet <- 1

for (p in 1:length(pthreslist)){

pthres <- pthreslist[p]
print("------------------")
print(pthres)

dat <- keep 
dat <- dat[dat$pval < pthres & is.na(dat$pval) == FALSE,]
exposure <- as.character(dat[1,"Phenotype"])
print(exposure)


tryCatch({

# Create exposure data set
exp_dat <- format_data(dat, type="exposure")
print(nrow(exp_dat))

# Clumping
exp_dat_clump <- clump_data(exp_dat)
print(paste0("Number of SNPs in exposure file after clumping: ",nrow(exp_dat_clump)))

# Read output data (no LD proxy possible so far!!!)
outcome_dat <- read_outcome_data(snps=exp_dat_clump$SNP, filename=outcomefilename,sep="\t")
outcome <- outcome_dat[1,"outcome"]

}, error = function(e) {
})


tryCatch({

# Harmonise data
dat <- harmonise_data(exposure_dat = exp_dat_clump, outcome_dat = outcome_dat)

# Perform MR
res <- mr(dat)
res$pthresMR <- pthres
res$GenExpDS <- genexpds
res$JLIM <- jlim
if (firstres == 1){
outres <- res
firstres <- 0
}else{outres <- rbind(outres,res)}


# Single SNP analysis
res_single <- mr_singlesnp(dat)

# Scatter plot
p1 <- mr_scatter_plot(res, dat)
# Forest plot
p2 <- mr_forest_plot(res_single)

#png(paste0(outputfolder,"/",exposure,"to",outcome,"_MRpthres_",pthres,"_plots.png"),width=15,height=10,units="cm",res=300,type="cairo")
#grid.arrange(p1[[1]],p2[[1]],nrow=1)
#dev.off()

}, error = function(e) {
})

tryCatch({
#png(paste0(outputfolder,"/",exposure,"to",outcome,"_MRpthres_",pthres,"_plots.png"),width=15,height=10,units="cm",res=300,type="cairo")
#grid.arrange(p1[[1]],p2[[1]],nrow=1)
#dev.off()

print(res)

}, error = function(e) {
})

# Horizontal Pleiotropy test
tryCatch({
plei <- mr_pleiotropy_test(dat)
plei$pthresMR <- pthres
plei$GenExpDS <- genexpds
plei$JLIM <- jlim
if (firstplei == 1){
outplei <- plei
firstplei <- 0
}else{outplei <- rbind(outplei,plei)}
}, warning = function(w) {
}, error = function(e) {
}, finally = {})

# Heterogeneity test
tryCatch({
het <- mr_heterogeneity(dat, method_list=c("mr_egger_regression", "mr_ivw"))
het$pthresMR <- pthres
het$GenExpDS <- genexpds
het$JLIM <- jlim
if (firsthet == 1){
outhet <- het
firsthet <- 0
}else{outhet <- rbind(outhet,het)}
}, warning = function(w) {
}, error = function(e) {
}, finally = {})

}

write.table(outres,paste0(outputfolder,"/",exposure,"to",outcome,"_MRresults.txt"),sep="\t",row.names=F,quote=F)

tryCatch({
write.table(outplei,paste0(outputfolder,"/",exposure,"to",outcome,"_MRpleiotropy.txt"),sep="\t",row.names=F,quote=F)
}, warning = function(w) {
}, error = function(e) {
}, finally = {})
tryCatch({
write.table(outhet,paste0(outputfolder,"/",exposure,"to",outcome,"_MRheterogeneity.txt"),sep="\t",row.names=F,quote=F)
}, warning = function(w) {
}, error = function(e) {
}, finally = {})

# MR Steiger directionality test
#out <- directionality_test(dat)
#out



