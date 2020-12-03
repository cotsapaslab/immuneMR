#! /bin/bash


PR_DIR=$1
WORK_DIR=$2

# Combine JLIM results
cat ${PR_DIR}/*_Tcells/loc*/2nd/loc*/jlim.out.tsv > ${PR_DIR}/jlim.out.tsv.Tcells 
cat ${PR_DIR}/*_Tcells/loc*/2nd/loc*/jlim.cfg.tsv > ${PR_DIR}/jlim.cfg.tsv.Tcells
cat ${PR_DIR}/*_Neutrophils/loc*/2nd/loc*/jlim.out.tsv > ${PR_DIR}/jlim.out.tsv.Neutrophils
cat ${PR_DIR}/*_Neutrophils/loc*/2nd/loc*/jlim.cfg.tsv > ${PR_DIR}/jlim.cfg.tsv.Neutrophils
cat ${PR_DIR}/*_Monocytes/loc*/2nd/loc*/jlim.out.tsv > ${PR_DIR}/jlim.out.tsv.Monocytes
cat ${PR_DIR}/*_Monocytes/loc*/2nd/loc*/jlim.cfg.tsv > ${PR_DIR}/jlim.cfg.tsv.Monocytes

cd ${PR_DIR}
cat jlim.out.tsv.Tcells jlim.out.tsv.Neutrophils jlim.out.tsv.Monocytes > jlim.out.tsv 
cat jlim.cfg.tsv.Tcells jlim.cfg.tsv.Neutrophils jlim.cfg.tsv.Monocytes > jlim.cfg.tsv

paste -d "\t" ${PR_DIR}/jlim.cfg.tsv ${PR_DIR}/jlim.out.tsv > ${PR_DIR}/jlim.res.tsv
awk '{ if (NR == 1 || $1 != "maintrID") print $0 }' ${PR_DIR}/jlim.res.tsv > ${PR_DIR}/jlim.res.tsv_2
mv ${PR_DIR}/jlim.res.tsv_2 ${PR_DIR}/jlim.res.tsv

#Prepare output file
rm -f ${PR_DIR}/jlim.info.tsv
touch ${PR_DIR}/jlim.info.tsv
printf "ImmPhen\tCond1stname\tCond1st\tCellType\tLocus\tGene\tCond2ndname\tCond2nd\n" > ${PR_DIR}/jlim.info.tsv

readarray lines < ${PR_DIR}/jlim.res.tsv
numlines=${#lines[@]}

for (( i=1; i<$numlines; i++ ))
do
index=$i

line="${lines[$index]}"
line=$(echo $line|tr -d '\n')

# Prepare relevant parameters

                OFS=$IFS
                IFS=$' '
                set -- $line
                array=( $@ )
                maintr=${array[0]}
                chr=${array[1]}
                genename=${array[9]}
		info=${array[10]}
                info2=${array[12]}
		IFS=$OFS

		OFS=$IFS
                IFS=$'.'
		set -- $info
		array=( $@ )
		genename2=${array[0]}
		gene="${genename}.${genename2}"
		chr=${array[1]}
		startbp=${array[2]}
		endbp=${array[3]}
		locus="locus.${chr}.${startbp}.${endbp}"
		cond2nd=${array[4]}
		IFS=$OFS

		delimiter="/2nd/"
		s=$info2$delimiter
		while [[ $s ]]; do
		    array1+=( "${s%%"$delimiter"*}" );
		    s=${s#*"$delimiter"};
		done;
		info3=${array1[0]}
		unset array1

		delimiter="/locus."
		s=$info3$delimiter
		while [[ $s ]]; do
                    array1+=( "${s%%"$delimiter"*}" );
                    s=${s#*"$delimiter"};
                done;
		info3=${array1[1]}
		info4=${array1[0]}
		unset array1

		OFS=$IFS
                IFS=$'.'
		set -- $info3
		array=( $@ )
		cond1st=${array[3]}
		IFS=$OFS

		OFS=$IFS
                IFS=$'_'
		set -- $info4
		array=( $@ )
		num=${#array[@]}
		num=$((num-1))
		CT=${array[$num]}
	 	IFS=$OFS


# Get info on 1st condition
if [ "$cond1st" == "uncond" ]; then
condition1="NA"
else
	trait=$maintr
	if [[ $trait == *"binary"* ]]
	then
	trait1=$trait
	else
	trait1="${trait}_norm"
	fi
	cond1file="${WORK_DIR}/MIP/GWAS/regression/${trait1}.${chr}.${startbp}.${endbp}.conds.txt"
	if [ "$cond1st" == "cond1" ]; then c1=1; fi
	if [ "$cond1st" == "cond2" ]; then c1=2; fi
	if [ "$cond1st" == "cond3" ]; then c1=3; fi
	if [ "$cond1st" == "cond4" ]; then c1=4; fi
	if [ "$cond1st" == "cond5" ]; then c1=5; fi
	if [ "$cond1st" == "cond6" ]; then c1=6; fi
	AKT_DIR=${PR_DIR}
	awk -v i="$c1" '{ if (NR != i) print $0}' $cond1file > ${AKT_DIR}/conds1st.txt
	cat ${AKT_DIR}/conds1st.txt | tr '\n' '_' > ${AKT_DIR}/conds1st.txt_2
	readarray condinfo1 < ${AKT_DIR}/conds1st.txt_2
	condinfo1=${condinfo1[0]}
	condinfo1=$(echo $condinfo1|tr -d '\n')
	condinfo1=${condinfo1::-1}
	rm ${AKT_DIR}/conds1st.txt_2
	rm ${AKT_DIR}/conds1st.txt
	condition1=$condinfo1
	unset $condinfo1	
fi

# Get info on 2nd condition
if [ "$cond2nd" == "uncond" ]; then
	condition2="NA"
else

	if [ "$cond2nd" == "cond1" ]; then condfile="${PR_DIR}/${CT}_assoc/${gene}.conds.txt_1"; fi
	if [ "$cond2nd" == "cond2" ]; then condfile="${PR_DIR}/${CT}_assoc/${gene}.conds.txt_2"; fi
	if [ "$cond2nd" == "cond3" ]; then condfile="${PR_DIR}/${CT}_assoc/${gene}.conds.txt_3"; fi	
	AKT_DIR=${PR_DIR}
	cat $condfile | tr '\n' '_' > ${AKT_DIR}/conds2nd.txt
        readarray condinfo < ${AKT_DIR}/conds2nd.txt
        condinfo=${condinfo[0]}
        condinfo=$(echo $condinfo|tr -d '\n')
        if [ "$condinfo" != "" ]; then
	condinfo=${condinfo::-1}
	else
	condinfo="NA"
	fi

        rm ${AKT_DIR}/conds2nd.txt
	condition2=$condinfo
	unset condinfo

fi
	
# Write info to file
printf "${maintr}\t${cond1st}\t${condition1}\t${CT}\t${locus}\t${gene}\t${cond2nd}\t${condition2}\n" >> ${PR_DIR}/jlim.info.tsv

done

# Combine files
paste -d "\t" ${PR_DIR}/jlim.info.tsv ${PR_DIR}/jlim.res.tsv > ${PR_DIR}/JLIMresults.txt























	
