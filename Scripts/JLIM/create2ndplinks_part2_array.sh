#! /bin/bash



#SBATCH --job-name=create2plink
#SBATCH --output=/home/cg859/scratch60/Logs/create2ndplink_part2_%A_%a_log.txt
#SBATCH --partition=general
#SBATCH --time=15:00:00
#SBATCH --array=2-201

WORK_DIR=$1
JLIM_DIR=$2
DATA_DIR=$3
phenoslist=$4
pThreshold2ndtrait=$5
pThreshold2ndtraitnum=$6
nperm=$7
numperm=$8



phenoslist1=$phenoslist

if [[ $phenoslist == "GenExp" ]]
then
	echo "GenExp phenos"
else
readarray phenos < $phenoslist
numphenos=${#phenos[@]}
fi


index=$((SLURM_ARRAY_TASK_ID-1))
#index=1
readarray lines < "${JLIM_DIR}/MS/MS-indexSNP.tsv"
line="${lines[$index]}"
line=$(echo $line|tr -d '\n')

                filename=${line}
                OFS=$IFS
                IFS=$' '
                set -- $filename
                array=( $@ )
                chr=${array[0]}
                startbp=${array[3]}
                endbp=${array[4]}
                IFS=$OFS

                if [[ $chr == "CHR" ]]
                then
                continue
                fi

		if  [[ $phenoslist1 == "GenExp" ]]
		then
		phenoslist="${DATA_DIR}/PhenoList${chr}.${startbp}.${endbp}.txt"
		readarray phenos < $phenoslist
		numphenos=${#phenos[@]}
		fi
			
		locus="locus.${chr}.${startbp}.${endbp}"
                echo $locus
 		
                plink --bfile ${JLIM_DIR}/2nd/${locus}/2nd --maf 0.05 --geno 0 --recode tab --out ${JLIM_DIR}/2nd/${locus}/2nd >/dev/null                
		

		for (( i=0; i<$numphenos; i++ ))
                do
                name="${phenos[$i]}"
	#	echo $name
		if [[ $name == *"binary"* ]]
		then
		bash ${WORK_DIR}/Scripts/JLIM/conditionalanalysis.sh ${WORK_DIR} ${JLIM_DIR}/2nd ${DATA_DIR} ${pThreshold2ndtrait} /${locus}/2nd ${name} logistic ${nperm} ${numperm}
		else
		bash ${WORK_DIR}/Scripts/JLIM/conditionalanalysis.sh ${WORK_DIR} ${JLIM_DIR}/2nd ${DATA_DIR} ${pThreshold2ndtrait} /${locus}/2nd ${name} linear ${nperm} ${numperm}
		fi
                done

