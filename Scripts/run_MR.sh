#! /bin/bash

#SBATCH --job-name=MRrun
#SBATCH --output=/home/cg859/scratch60/Logs/run_MR_%A_%a_log.txt
#SBATCH --partition=scavenge
#SBATCH --time=05:00:00
#SBATCH --array=1-16653%100

module load R/3.5.0-foss-2016b-avx2

liste=$1
PR_DIR=$2
WORK_DIR=$3

readarray MRlist < $liste
numMRlist=${#MRlist[@]}


i=$SLURM_ARRAY_TASK_ID
#i=15429

#for (( i=0; i<$numMRlist; i++ ))
#do

info="${MRlist[$i]}"
info=$(echo $info|tr -d '\n')

OFS=$IFS
IFS=$' '
set -- $info
array=( $@ )
GenExpDS=${array[1]}
Gene=${array[2]}
ImmPhen=${array[0]}
JLIM=${array[12]}
IFS=$OFS

if [[ $ImmPhen != *"binary"* ]] 
then
ImmPhen="${ImmPhen}_norm"
plinktest="linear"
else
plinktest="logistic"
fi

# Run MR analysis

if [ -f ${PR_DIR}/../MR/TwoSampleMR/${GenExpDS}_${Gene}to${ImmPhen}_MRresults.txt ]
then
printf ""
else
echo "$i---${GenExpDS}---${Gene}---${ImmPhen}---${JLIM}---"
Rscript ${WORK_DIR}/Scripts/TwoSampleMR.R ${PR_DIR}/${GenExpDS}_${Gene}.assoc.linear_prepared ${WORK_DIR}/MIP/GWAS/regression/${ImmPhen}/${ImmPhen}.${plinktest}.ADD_prepared ${GenExpDS} ${JLIM} ${PR_DIR}/../MR/TwoSampleMR
fi

#done


