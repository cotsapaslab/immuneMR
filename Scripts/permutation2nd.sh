#! /bin/bash

#SBATCH --partition=general
#SBATCH --time=20:00:00
#SBATCH --job-name=PRS
#SBATCH --output=/home/cg859/scratch60/Logs/PRS_permJLIM_%A_%a_log
#SBATCH --array=2583,2671,2677,2769,2770,2774,2801,2861,2860
#0-2864


PR_DIR=$1
plinkfile=$2
WORK_DIR=$3

module load PLINK/1.90-beta5.3
module load R

readarray lines < ${PR_DIR}/list_1stand2ndtraits.txt_unique2nd
numlines=${#lines[@]}
numperm=100000
nperm=$((numperm+100))

index=$((SLURM_ARRAY_TASK_ID*10))
endindex=$((index+10))

for (( p=$index; p<$endindex; p++ ))
do

line="${lines[$p]}"
line=$(echo $line|tr -d '\n')
OFS=$IFS
IFS=$' '
set -- $line
array=( $@ )
CT=${array[0]}
trait=${array[1]}
cond1st=${array[2]}
locus=${array[3]}
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


cp ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/${gene}.${chr}.${startbp}.${endbp}.${cond2nd}.assoc.linear.gz ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/${gene}.${chr}.${startbp}.${endbp}.${cond2nd}.assoc.linear_backup.gz
gunzip ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/${gene}.${chr}.${startbp}.${endbp}.${cond2nd}.assoc.linear.gz
awk '{ print $2 }' ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/${gene}.${chr}.${startbp}.${endbp}.${cond2nd}.assoc.linear > ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/${gene}.${chr}.${startbp}.${endbp}.${cond2nd}.snps.txt

if [ "$cond2nd" != "uncond" ];then
cat ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/${gene}.${chr}.${startbp}.${endbp}.${cond2nd}.snps.txt $condfile > ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/${gene}.${chr}.${startbp}.${endbp}.${cond2nd}.snps.txt_2
mv ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/${gene}.${chr}.${startbp}.${endbp}.${cond2nd}.snps.txt_2 ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/${gene}.${chr}.${startbp}.${endbp}.${cond2nd}.snps.txt
fi
 
plink --bfile $plinkfile --extract ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/${gene}.${chr}.${startbp}.${endbp}.${cond2nd}.snps.txt --make-bed --out ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/2nd
plink --bfile ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/2nd --recode tab --out ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/2nd

if [ -f ${PR_DIR}/${CT}_assoc/${gene}.${chr}.${startbp}.${endbp}.${cond2nd}.mperm.dump.all.gz ]; then
	echo "FOUND PERMUTATION FILE"
	cp ${PR_DIR}/${CT}_assoc/${gene}.${chr}.${startbp}.${endbp}.${cond2nd}.mperm.dump.all.gz ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/
	if [ "$cond2nd" == "uncond" ]; then
        plink --bfile ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/2nd --linear sex --pheno $phenofile --pheno-name $gene --covar $covarfile --ci 0.95 --out ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/${gene}.${chr}.${startbp}.${endbp}.${cond2nd}
        else
        plink --bfile ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/2nd --linear sex --pheno $phenofile --pheno-name $gene --covar $covarfile --ci 0.95 --condition-list $condfile --out ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/${gene}.${chr}.${startbp}.${endbp}.${cond2nd}
        fi
	awk '{ if (NR == 1 || $5 == "ADD") print $0 }' ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/${gene}.${chr}.${startbp}.${endbp}.${cond2nd}.assoc.linear > ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/${gene}.${chr}.${startbp}.${endbp}.${cond2nd}.assoc.linear.ADD
else
	echo "CREATE NEW PERMUTATION FILE"
	if [ "$cond2nd" == "uncond" ]; then
	plink --bfile ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/2nd --linear sex mperm=${nperm} --mperm-save-all --pheno $phenofile --pheno-name $gene --covar $covarfile --ci 0.95 --out ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/${gene}.${chr}.${startbp}.${endbp}.${cond2nd}
	else
	plink --bfile ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/2nd --linear sex mperm=${nperm} --mperm-save-all --pheno $phenofile --pheno-name $gene --covar $covarfile --ci 0.95 --condition-list $condfile --out ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/${gene}.${chr}.${startbp}.${endbp}.${cond2nd}
	fi
	awk '{ if (NR == 1 || $5 == "ADD") print $0 }' ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/${gene}.${chr}.${startbp}.${endbp}.${cond2nd}.assoc.linear > ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/${gene}.${chr}.${startbp}.${endbp}.${cond2nd}.assoc.linear.ADD
	Rscript ${WORK_DIR}/Scripts/checkperm.R ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/${gene}.${chr}.${startbp}.${endbp}.${cond2nd}.mperm.dump.all ${numperm} ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/${gene}.${chr}.${startbp}.${endbp}.${cond2nd}.assoc.linear.ADD
	gzip ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/${gene}.${chr}.${startbp}.${endbp}.${cond2nd}.mperm.dump.all
	cp ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/${gene}.${chr}.${startbp}.${endbp}.${cond2nd}.mperm.dump.all.gz ${PR_DIR}/${CT}_assoc/
	rm ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/${gene}.${chr}.${startbp}.${endbp}.${cond2nd}.assoc.linear.mperm	
	rm ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/${gene}.${chr}.${startbp}.${endbp}.${cond2nd}.mperm.dump.all.gz
fi

awk '{ if (NR == 1 || $12 != "NA") print $0 }' ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/${gene}.${chr}.${startbp}.${endbp}.${cond2nd}.assoc.linear.ADD > ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/${gene}.${chr}.${startbp}.${endbp}.${cond2nd}.assoc.linear

awk '{if (NR != 1) print $2 }' ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/${gene}.${chr}.${startbp}.${endbp}.${cond2nd}.assoc.linear > ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/${gene}.${chr}.${startbp}.${endbp}.${cond2nd}.snps.txt
plink --bfile $plinkfile --extract ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/${gene}.${chr}.${startbp}.${endbp}.${cond2nd}.snps.txt --make-bed --out ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/2nd
plink --bfile ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/2nd --recode tab --out ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/2nd

gzip ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/${gene}.${chr}.${startbp}.${endbp}.${cond2nd}.assoc.linear
gzip ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/2nd.ped
rm ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/2nd.b*
rm ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/2nd.fam
rm ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/2nd.log
rm ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/${gene}.${chr}.${startbp}.${endbp}.${cond2nd}.snps.txt
rm ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/${gene}.${chr}.${startbp}.${endbp}.${cond2nd}.log

OFS=$IFS
IFS=$'.'
set -- $gene
array=( $@ )
genename=${array[0]}
IFS=$OFS
mv ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/2nd.ped.gz ${PR_DIR}/${trait}_${CT}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${gene}.${cond2nd}/${genename}.ped.gz
IFS=$OFS


done


