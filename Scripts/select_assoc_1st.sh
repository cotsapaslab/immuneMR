#!/bin/bash

PR_DIR=$1

pthres=0.00001

rm -f ${PR_DIR}/list_1sttraits.txt
touch ${PR_DIR}/list_1sttraits.txt


readarray lines < ${PR_DIR}/list.txt
numlines=${#lines[@]}

for (( p=1; p<$numlines; p++ ))
do

line="${lines[$p]}"
line=$(echo $line|tr -d '\n')
OFS=$IFS
IFS=$' '
set -- $line
array=( $@ )
trait=${array[0]}
CT=${array[1]}
IFS=$OFS

ls ${PR_DIR}/${trait}_${CT}/1st/*.txt > ${PR_DIR}/${trait}_${CT}/assocs1st.txt
awk '{ if ($1 ~ /cond/) print $0 }' ${PR_DIR}/${trait}_${CT}/assocs1st.txt > ${PR_DIR}/${trait}_${CT}/assocs1st.txt_2
mv ${PR_DIR}/${trait}_${CT}/assocs1st.txt_2 ${PR_DIR}/${trait}_${CT}/assocs1st.txt

readarray conds < ${PR_DIR}/${trait}_${CT}/assocs1st.txt
numconds=${#conds[@]}


for (( i=0; i<$numconds; i++ ))
do


	file="${conds[$i]}"
	file=$(echo $file|tr -d '\n')
	filename=$file

	delimiter="1st"
	s=$file$delimiter
	while [[ $s ]]; do
	array1+=( "${s%%"$delimiter"*}" );
	s=${s#*"$delimiter"};
	done;
	file=${array1[1]}
	unset array1

	OFS=$IFS
        IFS=$'.'
        set -- $file
        array=( $@ )
        chr=${array[2]}
        startbp=${array[3]}
        endbp=${array[4]}
	cond=${array[5]}
	trait="${array[0]}.${array[1]}"
        IFS=$OFS
	
	# Get minimal p-value (only create JLIM folder if p < 1x10-5
	LC_ALL=C sort -k12 -g $filename > ${filename}_sort
	lead=`awk -v p=$pthres '{ if(NR == 2 && $12<p) print $2}' ${filename}_sort`
	rm ${filename}_sort
	if [ "$lead" != "" ]
	then
	mkdir ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond}
	mkdir ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond}/1st
	mkdir ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond}/2nd
 	cp $filename ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond}/1st/${trait}.${chr}.${startbp}.${endbp}.txt	
	awk -v chr="$chr" -v sbp="$startbp" -v ebp="$endbp" '{ if (NR == 1 || ($1 == chr && $4 == sbp && $5 == ebp)) print $0 }' ${PR_DIR}/${trait}_${CT}/1st/${trait}-index.tsv > ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond}/1st/${trait}-index.tsv
	printf "${CT}\t${trait}\tlocus${chr}.${startbp}.${endbp}\t${cond}\n" >> ${PR_DIR}/list_1sttraits.txt
	fi
	unset $lead



done

done


