args <- commandArgs(TRUE)

mat <- read.table(args[1],sep="\t",h=T)
snps <- read.table(args[2],sep="")

mat[mat$GENES %in% snps == TRUE,"BETA_GWAS"] <- mat[mat$GENES %in% snps == TRUE,"BETA_GWAS"] *-1

print(nrow(mat[mat$GENES %in% snps$V1 == TRUE,]))
if (nrow(mat[mat$GENES %in% snps$V1 == TRUE,]) > 0){
print("!!!!!switched")
}

write.table(mat,args[1],sep="\t",row.names=F,quote=F)

