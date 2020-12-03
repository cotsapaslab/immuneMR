#! /bin/bash

#SBATCH --partition=general
#SBATCH --time=1:00:00
#SBATCH --job-name=PRS
#SBATCH --output=/home/cg859/scratch60/Logs/PRS_run_JLIM_%A_%a_log
#SBATCH --array=2-494

module load PLINK/1.90-beta4.6
module load R
module load Pandoc
module load VCFtools
module load tabix

WORK_DIR=$2
PR_DIR=$3
SCR_DIR=$4

readarray lines < $1
numlines=${#lines[@]}

index=$((SLURM_ARRAY_TASK_ID-1))
#index=1


line="${lines[$index]}"
line=$(echo $line|tr -d '\n')
OFS=$IFS
IFS=$' '
set -- $line
array=( $@ )

firstname=${array[0]}
secondname=${array[1]}
assoc1dir=${array[2]}
if [[ $firstname == *"binary"* ]]
then
assoc1file="${assoc1dir}/${firstname}/${firstname}_results.RDS"
topassoc1file="${assoc1dir}/${firstname}/${firstname}_topresults_unique_1-5.txt"
else
assoc1file="${assoc1dir}/${firstname}_norm/${firstname}_norm_results.RDS"
topassoc1file="${assoc1dir}/${firstname}_norm/${firstname}_norm_topresults_unique_1-5.txt"
fi
secondphenfile=${array[3]}
secondphenos=${array[4]}
genotypefile=${array[5]}
DATA_DIR=${array[6]}
IFS=$OFS


#mkdir ${PR_DIR}/${firstname}_${secondname} 2>/dev/null
#mkdir ${PR_DIR}/${firstname}_${secondname}/1st 2>/dev/null
#mkdir ${PR_DIR}/${firstname}_${secondname}/2nd 2>/dev/null
AKT_DIR="${PR_DIR}/${firstname}_${secondname}"


# wenn Skript nur teilweise aufgeführt wird - check ob Locus vorhanden
if [ -f ${AKT_DIR}/1st/${firstname}-index.tsv ]
then
echo "ok"
else
exit
fi

# get list of genes to test
echo "Rscript ${SCR_DIR}/genesperlocus.R ${AKT_DIR}/1st/${firstname}-index.tsv ${secondphenfile} ${AKT_DIR}/2nd/"
Rscript ${SCR_DIR}/genesperlocus.R ${AKT_DIR}/1st/${firstname}-index.tsv ${secondphenfile} ${AKT_DIR}/2nd/


# Create reference LD file
${WORK_DIR}/Scripts/JLIM/jlim-master/bin/fetch.refld0.EUR.pl /home/cg859/project/BU_Deletion/JLIM_reference ${AKT_DIR}/1st/${firstname}-index.tsv ${WORK_DIR}/PRS-Project/ld0/

