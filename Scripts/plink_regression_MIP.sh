#! /usr/bin/bash

#SBATCH --job-name=GWAS_regression
#SBATCH --output=/home/cg859/scratch60/Logs/GWAS_regression_%A_%a_log.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00
#SBATCH --array=1-22
#SBATCH --partition=general


module load PLINK/1.90-beta5.3

dir=$1
plinkfile=$2
data_dir=$3
phenoslist=$4
chr=$SLURM_ARRAY_TASK_ID

readarray phenos < $phenoslist
numphenos=${#phenos[@]}

for (( i=0; i<$numphenos; i++ ))
do
name="${phenos[$i]}"
name=${name//[$'\t\r\n']}
mkdir ${dir}/regression/${name}
mkdir ${dir}/regression/${name}/chromosomes

	if [[ $name == *"binary"* ]]
        then
		plink --bfile ${plinkfile} --logistic sex --pheno ${data_dir}/phenotypes_transformed.txt --pheno-name ${name} --chr ${chr} --covar ${data_dir}/covariates_selected.txt --out ${dir}/regression/${name}/chromosomes/chr${chr} 
	else
		plink --bfile ${plinkfile} --linear sex --pheno ${data_dir}/phenotypes_transformed.txt --pheno-name ${name} --chr ${chr} --covar ${data_dir}/covariates_selected.txt --out ${dir}/regression/${name}/chromosomes/chr${chr}
	fi

done


