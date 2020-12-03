#! /bin/bash


module load PLINK/1.90-beta4.6

assocfile=$1
pthres=$2
genfile=$3
pheno=$4
phenoname=$5
covs=$6
listSNPs=$7


rm -f $listSNPs
touch $listSNPs

awk '{ if (NR > 1) print $2 }' ${assocfile} > ${assocfile}.snps

LC_ALL=C sort -k9 -g ${assocfile} > ${assocfile}.sort
lead=`awk -v p=$pthres '{ if(NR == 2 && $9<p) print $2}' ${assocfile}.sort`
awk '{ if (NR == 2) print $1"\t"$3"\t"$2 }' ${assocfile}.sort > $listSNPs
printf "${lead}\n" > ${listSNPs}_cond

while [ "$lead" != "" ]
do

	plink --bfile ${genfile} --extract ${assocfile}.snps --linear sex --pheno $pheno --pheno-name $phenoname --covar $covs --condition-list ${listSNPs}_cond --out ${assocfile}.cond
	awk '{if (NR == 1 || ($5 == "ADD" && $9 != "NA") ) print $0}' ${assocfile}.cond.assoc.linear > ${assocfile}.cond.assoc.linear.ADD
	LC_ALL=C sort -k9 -g ${assocfile}.cond.assoc.linear.ADD > ${assocfile}.cond.assoc.linear.ADD.sort
	lead=`awk -v p=$pthres '{ if(NR == 2 && $9<p) print $2}' ${assocfile}.cond.assoc.linear.ADD.sort`
	if [ "$lead" != "" ]
	then
	awk '{ if (NR == 2) print $1"\t"$3"\t"$2 }' ${assocfile}.cond.assoc.linear.ADD.sort >> $listSNPs
	printf "${lead}\n" >> ${listSNPs}_cond
	fi
	rm ${assocfile}.cond.*
done


rm ${listSNPs}_cond





