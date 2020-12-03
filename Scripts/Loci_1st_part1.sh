#! /bin/bash


#SBATCH --partition=general
#SBATCH --time=3:00:00
#SBATCH --job-name=PRS
#SBATCH --output=/home/cg859/scratch60/Logs/PRS_run_JLIM_%A_%a_log
#SBATCH --array=2-493

module load R

WORK_DIR=$2
PR_DIR=$3
SCR_DIR=$4

readarray lines < $1
numlines=${#lines[@]}

index=$((SLURM_ARRAY_TASK_ID-1))

line="${lines[$index]}"
line=$(echo $line|tr -d '\n')
echo $line
OFS=$IFS
IFS=$' '
set -- $line
array=( $@ )

firstname=${array[0]}
secondname=${array[1]}
assoc1dir=${array[2]}
if [[ $firstname == *"binary"* ]]
then
assoc1file="${assoc1dir}/${firstname}/${firstname}.logistic"
topassoc1file="${assoc1dir}/${firstname}/${firstname}.leads"
else
assoc1file="${assoc1dir}/${firstname}_norm/${firstname}_norm.linear"
topassoc1file="${assoc1dir}/${firstname}_norm/${firstname}_norm.leads"
fi
secondphenfile=${array[3]}
secondphenos=${array[4]}
genotypefile=${array[5]}
DATA_DIR=${array[6]}
IFS=$OFS


mkdir ${PR_DIR}/${firstname}_${secondname} 2>/dev/null
mkdir ${PR_DIR}/${firstname}_${secondname}/1st 2>/dev/null
mkdir ${PR_DIR}/${firstname}_${secondname}/2nd 2>/dev/null
AKT_DIR="${PR_DIR}/${firstname}_${secondname}"


# get summary statistics and index file for 1st trait (with just one GWAS hit)
awk '{ if (NR == 1 || $5 == "ADD") print $0}' ${assoc1file} > ${assoc1file}.ADD.${index}
Rscript ${SCR_DIR}/getassoc1.R ${topassoc1file} ${assoc1file}.ADD.${index} ${firstname} ${AKT_DIR}
rm ${assoc1file}.ADD.${index}

# MHC associations are filtered out in getassoc1.R, remove Files for those
if [ -f ${AKT_DIR}/${firstname}-index.tsv ]
then
mv ${AKT_DIR}/*.txt ${AKT_DIR}/1st/
mv ${AKT_DIR}/*.tsv ${AKT_DIR}/1st/
else
rm -r ${AKT_DIR}
exit
fi

echo "${AKT_DIR}/1st/${firstname}-index.tsv"

