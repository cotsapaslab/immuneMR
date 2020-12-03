#! /usr/bin/bash

module load R
module load PLINK/1.90-beta5.3

plinkfile=$1
phenofile=$2
covarfile=$3
pheno=$4
chr=$5
startbp=$6
endbp=$7
dir=$8
WD=$9

pthres=0.001
i=0

rm -f  ${dir}/${pheno}.conds.txt
touch ${dir}/${pheno}.conds.txt
plinktest="linear"

# Run first unconditional analysis to get first lead
plink --bfile ${plinkfile} --linear sex --pheno ${phenofile} --pheno-name ${pheno} --covar ${covarfile} --chr $chr --from-bp $startbp --to-bp $endbp --out ${dir}/${pheno}.${i}

grep ADD ${dir}/${pheno}.${i}.assoc.${plinktest} > ${dir}/${pheno}.${i}.assoc.${plinktest}.ADD
LC_ALL=C sort -k9 -g ${dir}/${pheno}.${i}.assoc.${plinktest}.ADD > ${dir}/${pheno}.${i}.assoc.${plinktest}.ADD.sort
lead=`awk -v p=$pthres '{ if(NR == 1 && $9<p) print $2}' ${dir}/${pheno}.${i}.assoc.${plinktest}.ADD.sort`

plink --bfile ${plinkfile} --chr $chr --from-bp $startbp --to-bp $endbp --make-bed --out ${dir}/${pheno}.subset
cp ${dir}/${pheno}.subset.bim ${dir}/${pheno}.snps.txt

while [ "$lead" != "" ] && [ $i -lt 3 ]
do
i=$((i+1))

printf "${lead}\n" >> ${dir}/${pheno}.conds.txt
plink --bfile ${plinkfile} --${plinktest} sex --pheno ${phenofile} --pheno-name ${pheno} --covar ${covarfile} --chr $chr --from-bp $startbp --to-bp $endbp --condition-list ${dir}/${pheno}.conds.txt --out ${dir}/${pheno}.${i}
grep ADD ${dir}/${pheno}.${i}.assoc.${plinktest} > ${dir}/${pheno}.${i}.assoc.${plinktest}.ADD
LC_ALL=C sort -k9 -g ${dir}/${pheno}.${i}.assoc.${plinktest}.ADD > ${dir}/${pheno}.${i}.assoc.${plinktest}.ADD.sort
grep -v NA ${dir}/${pheno}.${i}.assoc.${plinktest}.ADD.sort > ${dir}/${pheno}.${i}.assoc.${plinktest}.ADD.sort.nNA
# create LD file (in phase information) for lead SNP
plink --bfile ${dir}/${pheno}.subset --r2 in-phase --ld-snp $lead --out ${dir}/${pheno}.snps.r2
# Only select SNPs not in LD with all the leads
Rscript ${WD}/Scripts/selectSNPsbyR2.R ${dir}/${pheno}.snps.txt ${dir}/${pheno}.snps.r2.ld ${dir}/${pheno}.${i}.assoc.${plinktest}.ADD.sort.nNA
lead=`awk -v p=$pthres '{ if(NR == 1 && $9<p) print $2}' ${dir}/${pheno}.${i}.assoc.${plinktest}.ADD.sort.nNA`
echo "${i} --- ${lead} ---- !!!!!!!!!!!"
done

rm ${dir}/${pheno}.*.log
rm ${dir}/${pheno}.*.assoc*
rm ${dir}/${pheno}.subset.*
rm ${dir}/${pheno}.snps.txt
rm ${dir}/${pheno}.snps.r2.ld
