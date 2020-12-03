#! /bin/bash

#SBATCH --partition=general
#SBATCH --time=15:00:00 

PR_DIR=$1
typea=$2

if [[ $typea == *"woJLIM"* ]]
then
ls ${PR_DIR}/*2pheno_woJLIMlocus.prsice > ${PR_DIR}/listFiles.txt
else
	if [[ $typea == *"condJLIM"* ]]
	then
	ls ${PR_DIR}/*2pheno_condJLIMlocus.prsice > ${PR_DIR}/listFiles.txt
	else
	ls ${PR_DIR}/*2pheno_all.prsice > ${PR_DIR}/listFiles.txt
	fi
fi

readarray files < ${PR_DIR}/listFiles.txt
numfiles=${#files[@]}
echo $numfiles

for pthres in "5e-08" "1e-07" "5e-07" "1e-06" "5e-06" "1e-05" "5e-05" "0.0001" "0.001" "0.01"
do

echo $pthres

rm -f ${PR_DIR}/PRSresults.txt
touch ${PR_DIR}/PRSresults.txt

printf "DataSet\tGene\tImmPhen\tlocus\tSet\tThreshold\tPRS.R2\tP\tCoefficient\tStandard.Error\tNum_SNP\n" >> ${PR_DIR}/PRSresults.txt

for ((i = 0; i<$numfiles; i++))
do

filename="${files[$i]}"
OFS=$IFS
IFS=$'/'
set -- $filename
array=( $@ )
numarray=${#array[@]}
numarray=$((numarray-1))
file=${array[$numarray]}
file=$(echo $file|tr -d '\n')
#echo $file

IFS="\\_"
set -- $file
array=( $@ )
dataset=${array[0]}
phenotype1=${array[1]}
numarray=${#array[@]}
endpheno2=$((numarray-3))
#echo $file
	phenotype2=""
	for ((j = 2; j<$endpheno2; j++))
	do
	if  [ $j = 2 ]
	then
	phenotype2="${array[$j]}"
	else
	phenotype2=${phenotype2}"_"${array[$j]}
	fi
	done
locus="locus.${array[$endpheno2]}"
#echo "${dataset}---${phenotype1}---${phenotype2}---${locus}"
IFS=$OFS
awk -v th="$pthres" '{ if ($2 == th) print $0}' ${PR_DIR}/$file > ${PR_DIR}/${file}_2 
num_of_lines=$(wc -l <"${PR_DIR}/${file}_2")
if [[ $num_of_lines != "0" ]]
then
readarray infos < ${PR_DIR}/${file}_2
printf "${dataset}\t${phenotype1}\t${phenotype2}\t${locus}\t${infos[0]}" >> ${PR_DIR}/PRSresults.txt
fi
rm ${PR_DIR}/${file}_2
done

mv ${PR_DIR}/PRSresults.txt ${PR_DIR}/PRSresults_${typea}_${pthres}.txt


done
