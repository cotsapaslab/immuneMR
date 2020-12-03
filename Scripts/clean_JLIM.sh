#! /bin/bash

#SBATCH --partition=general
#SBATCH --time=12:00:00
#SBATCH --job-name=PRS
#SBATCH --output=/home/cg859/scratch60/Logs/PRS_permJLIM_%A_%a_log
#SBATCH --array=1-174,181-360


PR_DIR=$1
WORK_DIR=$2

module load PLINK/1.90-beta5.3
module load R

readarray lines < ${PR_DIR}/list_1stand2ndtraits.txt
numlines=${#lines[@]}


for (( p=0; p<$numlines; p++ ))
do

line="${lines[$p]}"
line=$(echo $line|tr -d '\n')
OFS=$IFS
IFS=$' '
set -- $line
array=( $@ )
CT=${array[0]}
trait=${array[1]}
locus=${array[3]}
cond1st=${array[2]}
gene=${array[4]}
cond2nd=${array[5]}
trait=${trait/\//}
IFS=$OFS

OFS=$IFS
IFS=$'.'
set -- $locus
array=( $@ )
chr=${array[0]}
chr=${chr/locus/}
startbp=${array[1]}
endbp=${array[2]}
IFS=$OFS

if [ "$CT" == "Tcells" ]
then
phenofile="${WORK_DIR}/BP/Data/Pheno/phenotypes_transformed.txt"
covarfile="${WORK_DIR}/BP/Data/Pheno/covariates_selected.txt"
fi
if [ "$CT" == "Monocytes" ]
then
phenofile="${WORK_DIR}/BP_Mono/Data/Pheno/phenotypes_transformed.txt"
covarfile="${WORK_DIR}/BP_Mono/Data/Pheno/covariates_selected.txt"
fi
if [ "$CT" == "Neutrophils" ]
then
phenofile="${WORK_DIR}/BP_Neutro/Data/Pheno/phenotypes_transformed.txt"
covarfile="${WORK_DIR}/BP_Neutro/Data/Pheno/covariates_selected.txt"
fi

if [ "$cond2nd" == "cond1" ]; then condfile="${PR_DIR}/${CT}_assoc/${gene}.conds.txt_1"; fi
if [ "$cond2nd" == "cond2" ]; then condfile="${PR_DIR}/${CT}_assoc/${gene}.conds.txt_2"; fi
if [ "$cond2nd" == "cond3" ]; then condfile="${PR_DIR}/${CT}_assoc/${gene}.conds.txt_3"; fi


AKT_DIR="${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}"

if [ -f ${AKT_DIR}/jlim.cfg.tsv ];then
	if [ -f ${AKT_DIR}/jlim.out.tsv ]; then
	echo "ok"
	else
	mv ${AKT_DIR}/jlim.cfg.tsv ${AKT_DIR}/jlim.cfg.tsv_notworked
	fi
fi


done


