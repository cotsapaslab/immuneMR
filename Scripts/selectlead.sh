#! /bin/bash


PR_DIR=$1

rm -f ${PR_DIR}/JLIMresults_selected_leads.txt
touch ${PR_DIR}/JLIMresults_selected_leads.txt
printf "DataSet\tGene\tChr\tidxSNP\tidxBP\tStartBP\tEndBP\tConds\tLead\n" >> ${PR_DIR}/JLIMresults_selected_leads.txt

rm -f ${PR_DIR}/JLIMresults_selected_complete.txt
touch ${PR_DIR}/JLIMresults_selected_complete.txt
printf "maintrID\tDataSet\tGene\tChr\tidxSNP\tidxBP\tStartBP\tEndBP\tConds\tLead\tSTAT\tminp\tJLIM\n" >> ${PR_DIR}/JLIMresults_selected_complete.txt


readarray lines < ${PR_DIR}/JLIMresults_selected.txt
numlines=${#lines[@]}



for (( i=1; i<$numlines; i++ ))
do

echo "${i}---${numlines}"

line="${lines[$i]}"
line=$(echo $line|tr -d '\n')

                filename=${line}
                OFS=$IFS
                IFS=$' '
                set -- $filename
                array=( $@ )
                ImmPhen=${array[0]}
                CellType=${array[1]}
		Gene=${array[2]}
		chrom=${array[3]}
		startbp=${array[4]}
                endbp=${array[5]}
		idxSNP=${array[6]}
		idxBP=${array[7]}
		stat=${array[8]}
		minp=${array[9]}
		cond=${array[10]}
		jlim=${array[11]}
		jlimfdr=${array[12]}
		cond1st=${array[13]}
		cond2nd=${array[14]}
                IFS=$OFS


assfile="${PR_DIR}/${CellType}_assoc/${Gene}.${cond2nd}.${chrom}.${startbp}.${endbp}.assoc.txt"
LC_ALL=C sort -k12 -g $assfile > ${assfile}_sorted
lead=`awk -v p=$pthres '{ if(NR == 2) print $2}' ${assfile}_sorted`
rm ${assfile}_sorted

printf "${CellType}\t${Gene}\t${chrom}\t${idxSNP}\t${idxBP}\t${startbp}\t${endbp}\t${cond}\t${lead}\n" >> ${PR_DIR}/JLIMresults_selected_leads.txt
printf "${ImmPhen}\t${CellType}\t${Gene}\t${chrom}\t${idxSNP}\t${idxBP}\t${startbp}\t${endbp}\t${cond}\t${lead}\t${stat}\t${minp}\t${jlimfdr}\n" >>  ${PR_DIR}/JLIMresults_selected_complete.txt
unset lead

done

