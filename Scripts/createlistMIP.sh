#! /bin/bash 


phenoslist=$1  #/home/cg859/scratch60/MIP/Data/Pheno/ImmPhen/PhenoList.txt
dir=$2 # /home/cg859/scratch60/MIP/GWAS
outfile=$3 # /home/cg859/scratch60/PRS-Project/MIP_BIP/list.txt


assoc1dir="/home/cg859/scratch60/MIP/GWAS/regression/"
SecondPhen1Tc="/home/cg859/scratch60/BP/Data/Pheno/Gtex/all.expression.bed"
SecondPhen1Mo="/home/cg859/scratch60/BP_Mono/Data/Pheno/Gtex/all.expression.bed"
SecondPhen1Ne="/home/cg859/scratch60/BP_Neutro/Data/Pheno/Gtex/all.expression.bed"
SecondPhenTc="/home/cg859/scratch60/BP/Data/Pheno/phenotypes_transformed.txt"
SecondPhenMo="/home/cg859/scratch60/BP_Mono/Data/Pheno/phenotypes_transformed.txt"
SecondPhenNe="/home/cg859/scratch60/BP_Neutro/Data/Pheno/phenotypes_transformed.txt"
GenoFile="/home/cg859/scratch60/BP/Data/Geno/checkQC/all_noMiss_noHet_noRel_noOutl"
DataDirTc="/home/cg859/scratch60/BP/Data/Pheno"
DataDirMo="/home/cg859/scratch60/BP_Mono/Data/Pheno"
DataDirNe="/home/cg859/scratch60/BP_Neutro/Data/Pheno"

readarray phenos < $phenoslist
numphenos=${#phenos[@]}

rm -f ${outfile}
touch ${outfile}
printf "Firstname\tSecondname\tAssoc1dir\tSecondPhen1\tSecondPhen\tGenoFile\tDataDir\n" >> ${outfile}

for (( i=0; i<$numphenos; i++ ))
do

	name="${phenos[$i]}"
	name=${name//[$'\t\r\n']}
	sigfile="${dir}/regression/${name}/${name}.leads"
	name=$(echo "${name/_norm/}")
	numb=$(sed -n '$=' $sigfile)
	if (($numb > 1)); then
	echo $numb
	printf "${name}\tTcells\t${assoc1dir}\t${SecondPhen1Tc}\t${SecondPhenTc}\t${GenoFile}\t${DataDirTc}\n" >> ${outfile}
	printf "${name}\tMonocytes\t${assoc1dir}\t${SecondPhen1Mo}\t${SecondPhenMo}\t${GenoFile}\t${DataDirMo}\n" >> ${outfile}
	printf "${name}\tNeutrophils\t${assoc1dir}\t${SecondPhen1Ne}\t${SecondPhenNe}\t${GenoFile}\t${DataDirNe}\n" >> ${outfile}

fi

done




