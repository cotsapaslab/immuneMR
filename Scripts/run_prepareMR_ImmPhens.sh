#! /bin/bash

#SBATCH --job-name=prepImmPhensMR
#SBATCH --output=${WORK_DIR}/Logs/MR_prepareImmPhens_run_%A_%a_log.txt
#SBATCH --array=0,2-163
#SBATCH --partition=general


module load R

liste=$1
PR_DIR=$2
WORK_DIR=$3

readarray MRlist < $liste


i=$SLURM_ARRAY_TASK_ID
#i=1

ImmPhen="${MRlist[$i]}"
ImmPhen=$(echo $ImmPhen|tr -d '\n')

if [[ $ImmPhen == *"binary"* ]]
then
plinktest="logistic"
else
plinktest="linear"
fi

# Prepare outcome file 
if [ -f ${WORK_DIR}/MIP/GWAS/regression/${ImmPhen}/${ImmPhen}.${plinktest}.ADD_prepared ]
then
echo "Output file prepared"
else
echo "Prepare outcome file..."
Rscript ${WORK_DIR}/Scripts/prepareData_ImmPhen.R ${WORK_DIR}/MIP/GWAS/regression/${ImmPhen}/${ImmPhen}.${plinktest}.ADD ${WORK_DIR}/MIP/Data/Genotypes/checkQC/LabExMI_imputation_816x5699237_noMiss_noHet_noRel_noOutl_HWE_MAF_QCd.bim ${WORK_DIR}/MIP/Data/Genotypes/checkQC/LabExMI_imputation_816x5699237_noMiss_noHet_noRel_noOutl_HWE_MAF_QCd.frq ${ImmPhen} ${WORK_DIR}/MIP/GWAS/regression/${ImmPhen}/${ImmPhen}.${plinktest}.ADD_prepared
fi


