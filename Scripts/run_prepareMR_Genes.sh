#! /bin/bash


#SBATCH --job-name=prepGenesMR
#SBATCH --output=${WORK_DIR}/Logs/MR_prepareGenes_run_%A_%a_log.txt
#SBATCH --array=0-1106%500
#SBATCH --time=12:00:00
#SBATCH --partition=general

module load R

liste=$1
PR_DIR=$2
WORK_DIR=$3

readarray MRlist < $liste



startindex=$((SLURM_ARRAY_TASK_ID * 10))
startindex=$((startindex + 1))
endindex=$((startindex + 10))

for (( cpi=$startindex; cpi<$endindex; cpi++ ))
do


i=$cpi

echo $i
echo "number"

line="${MRlist[$i]}"
line=$(echo $line|tr -d '\n')


filename=${line}
                OFS=$IFS
                IFS=$'_'
                set -- $filename
                array=( $@ )
                DS=${array[0]}
                Gene=${array[1]}
                IFS=$OFS




# Prepare exposure file
if [ -f ${PR_DIR}/PRS/${DS}_${Gene}.assoc.linear_prepared ]
then
echo "Exposure file prepared"
else
echo "Prepare exposure file..."
gunzip ${PR_DIR}/PRS/${DS}_${Gene}.assoc.linear.gz
Rscript ${WORK_DIR}/Scripts/prepareData_neu.R ${PR_DIR}/PRS/${DS}_${Gene}.assoc.linear ${WORK_DIR}/BP/Data/Geno/checkQC/all_noMiss_noHet_noRel_noOutl.bim ${WORK_DIR}/BP/Data/Geno/checkQC/all_noMiss_noHet_noRel_noOutl.frq ${DS}_${Gene} ${PR_DIR}/PRS/${DS}_${Gene}.assoc.linear_prepared
gzip ${PR_DIR}/PRS/${DS}_${Gene}.assoc.linear
fi

done

