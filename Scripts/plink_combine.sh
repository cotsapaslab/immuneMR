#! /usr/bin/bash

#SBATCH --job-name=GWAS_combine
#SBATCH --output=/home/cg859/scratch60/Logs/GWAS_combine_%A_%a_log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00
#SBATCH --array=0-163
#SBATCH --partition=general


module load PLINK/1.90-beta5.3

wd=$1
dir=$2
phenoslist=$3
index=$SLURM_ARRAY_TASK_ID

readarray phenos < $phenoslist
numphenos=${#phenos[@]}

name="${phenos[$index]}"
name=${name//[$'\t\r\n']}

if [[ $name == *"binary"* ]]
then
	cp ${wd}/Scripts/GWAS/combine_logistic.sh ${dir}/regression/${name}/chromosomes/
        cd ${dir}/regression/${name}/chromosomes/
        bash ${dir}/regression/${name}/chromosomes/combine_logistic.sh ${name}
        cd ${wd}
else
	cp ${wd}/Scripts/GWAS/combine.sh ${dir}/regression/${name}/chromosomes/
	cd ${dir}/regression/${name}/chromosomes/
	bash ${dir}/regression/${name}/chromosomes/combine.sh ${name}
	cd ${wd}
fi

