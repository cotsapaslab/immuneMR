#! /usr/bin/bash


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

rm -f  ${dir}/${pheno}.${chr}.${startbp}.${endbp}.conds.txt
touch ${dir}/${pheno}.${chr}.${startbp}.${endbp}.conds.txt


# Run first unconditional analysis to get first lead
if [[ $pheno == *"binary"* ]]
then
plink --bfile ${plinkfile} --logistic sex --pheno ${phenofile} --pheno-name ${pheno} --covar ${covarfile} --chr $chr --from-bp $startbp --to-bp $endbp --out ${dir}/${pheno}.${chr}.${startbp}.${endbp}.${i}
plinktest="logistic"
else
plink --bfile ${plinkfile} --linear sex --pheno ${phenofile} --pheno-name ${pheno} --covar ${covarfile} --chr $chr --from-bp $startbp --to-bp $endbp --out ${dir}/${pheno}.${chr}.${startbp}.${endbp}.${i}
plinktest="linear"
fi

grep ADD ${dir}/${pheno}.${chr}.${startbp}.${endbp}.${i}.assoc.${plinktest} > ${dir}/${pheno}.${chr}.${startbp}.${endbp}.${i}.assoc.${plinktest}.ADD
LC_ALL=C sort -k9 -g ${dir}/${pheno}.${chr}.${startbp}.${endbp}.${i}.assoc.${plinktest}.ADD > ${dir}/${pheno}.${chr}.${startbp}.${endbp}.${i}.assoc.${plinktest}.ADD.sort
lead=`awk -v p=$pthres '{ if(NR == 1 && $9<p) print $2}' ${dir}/${pheno}.${chr}.${startbp}.${endbp}.${i}.assoc.${plinktest}.ADD.sort`

plink --bfile ${plinkfile} --chr $chr --from-bp $startbp --to-bp $endbp --make-bed --out ${dir}/${pheno}.${chr}.${startbp}.${endbp}.subset
cp ${dir}/${pheno}.${chr}.${startbp}.${endbp}.subset.bim ${dir}/${pheno}.${chr}.${startbp}.${endbp}.snps.txt

while [ "$lead" != "" ]
do
i=$((i+1))

printf "${lead}\n" >> ${dir}/${pheno}.${chr}.${startbp}.${endbp}.conds.txt
plink --bfile ${plinkfile} --${plinktest} sex --pheno ${phenofile} --pheno-name ${pheno} --covar ${covarfile} --chr $chr --from-bp $startbp --to-bp $endbp --condition-list ${dir}/${pheno}.${chr}.${startbp}.${endbp}.conds.txt --out ${dir}/${pheno}.${chr}.${startbp}.${endbp}.${i}
grep ADD ${dir}/${pheno}.${chr}.${startbp}.${endbp}.${i}.assoc.${plinktest} > ${dir}/${pheno}.${chr}.${startbp}.${endbp}.${i}.assoc.${plinktest}.ADD
LC_ALL=C sort -k9 -g ${dir}/${pheno}.${chr}.${startbp}.${endbp}.${i}.assoc.${plinktest}.ADD > ${dir}/${pheno}.${chr}.${startbp}.${endbp}.${i}.assoc.${plinktest}.ADD.sort
grep -v NA ${dir}/${pheno}.${chr}.${startbp}.${endbp}.${i}.assoc.${plinktest}.ADD.sort > ${dir}/${pheno}.${chr}.${startbp}.${endbp}.${i}.assoc.${plinktest}.ADD.sort.nNA
# create LD file (in phase information) for lead SNP
plink --bfile ${dir}/${pheno}.${chr}.${startbp}.${endbp}.subset --r2 in-phase --ld-snp $lead --out ${dir}/${pheno}.${chr}.${startbp}.${endbp}.snps.r2
# Only select SNPs not in LD with all the leads
Rscript ${WD}/Scripts/selectSNPsbyR2.R ${dir}/${pheno}.${chr}.${startbp}.${endbp}.snps.txt ${dir}/${pheno}.${chr}.${startbp}.${endbp}.snps.r2.ld ${dir}/${pheno}.${chr}.${startbp}.${endbp}.${i}.assoc.${plinktest}.ADD.sort.nNA
lead=`awk -v p=$pthres '{ if(NR == 1 && $9<p) print $2}' ${dir}/${pheno}.${chr}.${startbp}.${endbp}.${i}.assoc.${plinktest}.ADD.sort.nNA`

done

rm ${dir}/${pheno}.${chr}.${startbp}.${endbp}.*.log
rm ${dir}/${pheno}.${chr}.${startbp}.${endbp}.*.assoc*
rm ${dir}/${pheno}.${chr}.${startbp}.${endbp}.subset.*
rm ${dir}/${pheno}.${chr}.${startbp}.${endbp}.snps.txt
rm ${dir}/${pheno}.${chr}.${startbp}.${endbp}.snps.r2.ld
