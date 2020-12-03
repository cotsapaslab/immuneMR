#! /bin/bash

#SBATCH --job-name=PRSice_assoc_cond
#SBATCH --output=/home/cg859/scratch60/Logs/PRS_assoc_cond_%A_%a_log.txt
#SBATCH --partition=general
#SBATCH --cpus-per-task=1
#SBATCH --time=8:00:00
#SBATCH --mem-per-cpu=5000
#SBATCH --array=1-13255%200


module load PLINK/1.90-beta4.6
module load R


WORK_DIR=$1
PR_DIR=$2

readarray lines < ${PR_DIR}/../JLIMresults_selected_leads_unique.txt

numlines=${#lines[@]}

i=$SLURM_ARRAY_TASK_ID

line="${lines[$i]}"
line=$(echo $line|tr -d '\n')

filename=${line}
                OFS=$IFS
                IFS=$' '
                set -- $filename
                array=( $@ )
                dataset1=${array[0]}
                phenotype1=${array[1]}
		lead=${array[2]}
                IFS=$OFS

	
        if [ "$dataset1" == "Tcells" ]
        then
        pheno1="${WORK_DIR}/BP/Data/Pheno/phenotypes_transformed.txt"
	cov1="${WORK_DIR}/BP/Data/Pheno/covariates_selected.txt"
	fi

	if [ "$dataset1" == "Monocytes" ]
        then
        pheno1="${WORK_DIR}/BP_Mono/Data/Pheno/phenotypes_transformed.txt"
        cov1="${WORK_DIR}/BP_Mono/Data/Pheno/covariates_selected.txt"
        fi

	if [ "$dataset1" == "Neutrophils" ]
        then
        pheno1="${WORK_DIR}/BP_Neutro/Data/Pheno/phenotypes_transformed.txt"
        cov1="${WORK_DIR}/BP_Neutro/Data/Pheno/covariates_selected.txt"
        fi
	

gen1="${WORK_DIR}/BP/Data/Geno/checkQC/all_noMiss_noHet_noRel_noOutl"

# Run GWAS 
#if [ -f ${PR_DIR}/${dataset1}_${phenotype1}_cond${lead}.assoc.linear.gz ]; then
#echo "already done"
#else

plink --bfile ${gen1} --linear hide-covar sex --pheno ${pheno1} --pheno-name ${phenotype1} --covar ${cov1} --condition ${lead} --out ${PR_DIR}/${dataset1}_${phenotype1}_cond${lead} 

grep ADD ${PR_DIR}/${dataset1}_${phenotype1}_cond${lead}.assoc.linear > ${PR_DIR}/${dataset1}_${phenotype1}_cond${lead}.assoc.linear.ADD
rm ${PR_DIR}/${dataset1}_${phenotype1}_cond${lead}.assoc.linear
grep -v NA ${PR_DIR}/${dataset1}_${phenotype1}_cond${lead}.assoc.linear.ADD > ${PR_DIR}/${dataset1}_${phenotype1}_cond${lead}.assoc.linear.ADD.NA
rm ${PR_DIR}/${dataset1}_${phenotype1}_cond${lead}.assoc.linear.ADD

Rscript ${WORK_DIR}/Scripts/addA2.R ${PR_DIR}/gen1.bim ${PR_DIR}/${dataset1}_${phenotype1}_cond${lead}.assoc.linear.ADD.NA ${PR_DIR}/${dataset1}_${phenotype1}_cond${lead}.assoc.linear
rm ${PR_DIR}/${dataset1}_${phenotype1}_cond${lead}.assoc.linear.ADD.NA
awk '{ if (NR == 1 || $9 < 0.01) print $0}' ${PR_DIR}/${dataset1}_${phenotype1}_cond${lead}.assoc.linear > ${PR_DIR}/${dataset1}_${phenotype1}_cond${lead}.assoc.linear_2
mv ${PR_DIR}/${dataset1}_${phenotype1}_cond${lead}.assoc.linear_2 ${PR_DIR}/${dataset1}_${phenotype1}_cond${lead}.assoc.linear 
gzip ${PR_DIR}/${dataset1}_${phenotype1}_cond${lead}.assoc.linear
echo "done"
###fi
