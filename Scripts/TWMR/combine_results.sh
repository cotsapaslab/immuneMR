#! /bin/bash


PR_DIR=$1

rm -f ${PR_DIR}/all.alpha

ls ${PR_DIR}/*.alpha > ${PR_DIR}/listFiles.txt

readarray files < ${PR_DIR}/listFiles.txt
numfiles=${#files[@]}

touch ${PR_DIR}/all.alpha

printf "DataSet\tImmPhen\tJLIM\tGene\talpha\tSE\tP\tNsnps\tNgene\n" >> ${PR_DIR}/all.alpha


for ((i = 0; i<$numfiles; i++))
do


filename="${files[$i]}"
OFS=$IFS

IFS=$'/'
set -- $filename
array=( $@ )
numarray=${#array[@]}
numarray=$((numarray-1))
filename=${array[$numarray]}


str=$filename
delimiter="_EN"
s=$str$delimiter
array=();
while [[ $s ]]; do
    array+=( "${s%%"$delimiter"*}" );
    s=${s#*"$delimiter"};
done;
dataset=${array[0]}
rest=${array[1]}


IFS="\\_"
set -- $rest
array=( $@ )
phenotype=""
numarray=${#array[@]}
	for ((j = 1; j<$numarray; j++))
	do
	if  [ $j = 1 ]
	then
	phenotype="${array[$j]}"
	else
	phenotype=${phenotype}"_"${array[$j]}
	fi
	done

str=$phenotype
delimiter="_JLIM"
s=$str$delimiter
array=();
while [[ $s ]]; do
    array+=( "${s%%"$delimiter"*}" );
    s=${s#*"$delimiter"};
done;
phenotype=${array[0]}
jlim="JLIM${array[1]}"

str=$jlim
delimiter=".alp"
s=$str$delimiter
array=();
while [[ $s ]]; do
    array+=( "${s%%"$delimiter"*}" );
    s=${s#*"$delimiter"};
done;
jlim=${array[0]}

IFS=$OFS

readarray infos < ${PR_DIR}/${filename}

if [ "${infos[1]}" != "" ]
then  
printf "${dataset}\t${phenotype}\t${jlim}\t${infos[1]}" >> ${PR_DIR}/all.alpha
fi

done

