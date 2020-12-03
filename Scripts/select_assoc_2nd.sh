#! /bin/bash

PR_DIR=$1
pthres=0.001

rm -f ${PR_DIR}/list_1stand2ndtraits.txt
touch ${PR_DIR}/list_1stand2ndtraits.txt

readarray lines < ${PR_DIR}/list_1sttraits.txt
numlines=${#lines[@]}

for (( p=0; p<$numlines; p++ ))
do
echo $p
line="${lines[$p]}"
line=$(echo $line|tr -d '\n')
OFS=$IFS
IFS=$' '
set -- $line
array=( $@ )
CT=${array[0]}
trait=${array[1]}
locus=${array[2]}
cond1st=${array[3]}
trait=${trait/\//}
IFS=$OFS

OFS=$IFS
IFS=$'.'
set -- $locus
array=( $@ )
chr=${array[0]}
chr=${chr/locus/}
startbp=${array[1]}
endbp=${array[2]}
IFS=$OFS


readarray genes < ${PR_DIR}/${trait}_${CT}/2nd/PhenoList${chr}.${startbp}.${endbp}.txt
numgenes=${#genes[@]}

for (( g=0; g<$numgenes; g++))
do

	gene="${genes[$g]}"
	gene=$(echo $gene|tr -d '\n')
	for cond in uncond cond1 cond2 cond3
	do
		if [ -f ${PR_DIR}/${CT}_assoc/${gene}.${cond}.assoc.txt ]
		then
		awk -v chr="$chr" -v sbp="$startbp" -v ebp="$endbp" '{ if ( NR == 1 || ($1 == chr && $3 >= sbp && $3 <= ebp)) print $0 }' ${PR_DIR}/${CT}_assoc/${gene}.${cond}.assoc.txt > ${PR_DIR}/${CT}_assoc/${gene}.${cond}.${chr}.${startbp}.${endbp}.assoc.txt
		LC_ALL=C sort -k12 -g ${PR_DIR}/${CT}_assoc/${gene}.${cond}.${chr}.${startbp}.${endbp}.assoc.txt > ${PR_DIR}/${CT}_assoc/${gene}.${cond}.${chr}.${startbp}.${endbp}.assoc.txt_sort
      		lead=`awk -v p=$pthres '{ if(NR == 2 && $12<p) print $2}' ${PR_DIR}/${CT}_assoc/${gene}.${cond}.${chr}.${startbp}.${endbp}.assoc.txt_sort`
		rm ${PR_DIR}/${CT}_assoc/${gene}.${cond}.${chr}.${startbp}.${endbp}.assoc.txt_sort
		if [ "$lead" != "" ]
        	then
		mkdir ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond}
	 	cp ${PR_DIR}/${CT}_assoc/${gene}.${cond}.${chr}.${startbp}.${endbp}.assoc.txt ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond}/${gene}.${chr}.${startbp}.${endbp}.${cond}.assoc.linear
		gzip ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond}/${gene}.${chr}.${startbp}.${endbp}.${cond}.assoc.linear	
		printf "${CT}\t${trait}\t${cond1st}\t${locus}\t${gene}\t${cond}\n" >> ${PR_DIR}/list_1stand2ndtraits.txt
		fi
		unset $lead
		fi 
	done

done


done


