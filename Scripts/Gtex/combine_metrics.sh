#! /bin/bash



idsFile=$1
dir=$2


readarray ids < $idsFile
numids=${#ids[@]}

rm ${dir}/all.metrics.tsv
touch ${dir}/all.metrics.tsv


for (( i=0; i<$numids; i++ ))
do

	sample_id="${ids[$i]}"
	sample_id="$(echo "$sample_id"|tr -d '\n')"


	if (( $i == 1 ))
	then
		head -n 1 ${dir}/alldone/${sample_id}.metrics.tsv > ${dir}/alldone/header.metrics
	fi

	sed -n '2p' < ${dir}/alldone/${sample_id}.metrics.tsv >> ${dir}/all.metrics.tsv

done
cat ${dir}/alldone/header.metrics ${dir}/all.metrics.tsv > ${dir}/all.metrics.tsv_2
mv ${dir}/all.metrics.tsv_2 ${dir}/all.metrics.tsv

rm ${dir}/alldone/header.metrics

