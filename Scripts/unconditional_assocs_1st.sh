#! /bin/bash

#SBATCH --partition=scavenge
#SBATCH --time=2:00:00
#SBATCH --job-name=cond_1st
#SBATCH --output=/home/cg859/scratch60/Logs/PRS_conditionals1st_%A_%a_log
#SBATCH --array=1-1432

module load PLINK/1.90-beta5.3

readarray lines < $1
plinkfile=$2
phenofile=$3
covarfile=$4
dir=$5
pr_dir=$6

index=$((SLURM_ARRAY_TASK_ID-1))

line="${lines[$index]}"
line=$(echo $line|tr -d '\n')

OFS=$IFS
IFS=$' '
set -- $line
array=( $@ )
IFS=$OFS


num=${array[0]}
info=${array[1]}

echo $num
echo $info

if [[ $info == *"binary"* ]]
then
delimiter="_binary"
else
delimiter="_norm"
fi
delimiter1="regression/"
delimiter2="."

s=$info$delimiter
while [[ $s ]]; do
    array1+=( "${s%%"$delimiter"*}" );
    s=${s#*"$delimiter"};
done;
pheno=${array1[0]}
rest=${array1[1]}

s=$pheno$delimiter1
while [[ $s ]]; do
    array2+=( "${s%%"$delimiter1"*}" );
    s=${s#*"$delimiter1"};
done;
pheno=${array2[1]}


s=$rest$delimiter2
while [[ $s ]]; do
    array3+=( "${s%%"$delimiter2"*}" );
    s=${s#*"$delimiter2"};
done;
chr=${array3[1]}
startbp=${array3[2]}
endbp=${array3[3]}


if [[ $info == *"binary"* ]]
then
phenoold="${pheno}_binary"
pheno="${pheno}_binary"
plinktest="logistic"
else
phenoold=$pheno
plinktest="linear"
pheno="${pheno}_norm"
fi


plink --bfile ${plinkfile} --${plinktest} sex --pheno ${phenofile} --pheno-name ${pheno} --covar ${covarfile} --ci 0.95 --chr $chr --from-bp $startbp --to-bp $endbp --out ${dir}/${pheno}.${chr}.${startbp}.${endbp}.uncond
grep ADD ${dir}/${pheno}.${chr}.${startbp}.${endbp}.uncond.assoc.${plinktest} > ${dir}/${pheno}.${chr}.${startbp}.${endbp}.uncond.assoc.${plinktest}.ADD
grep -v NA ${dir}/${pheno}.${chr}.${startbp}.${endbp}.uncond.assoc.${plinktest}.ADD > ${dir}/${pheno}.${chr}.${startbp}.${endbp}.uncond.assoc.${plinktest}.ADD.nNA
printf "CHR\tSNP\tBP\tA1\tTEST\tNMISS\tBETA\tSE\tLxy\tHxy\tSTAT\tP\n" > ${dir}/${phenoold}.${chr}.${startbp}.${endbp}.uncond.txt_head
cat ${dir}/${phenoold}.${chr}.${startbp}.${endbp}.uncond.txt_head ${dir}/${pheno}.${chr}.${startbp}.${endbp}.uncond.assoc.${plinktest}.ADD.nNA > ${dir}/${phenoold}.${chr}.${startbp}.${endbp}.uncond.assoc.txt
cp ${dir}/${phenoold}.${chr}.${startbp}.${endbp}.uncond.assoc.txt ${pr_dir}/${phenoold}_Tcells/1st/${phenoold}.${chr}.${startbp}.${endbp}.uncond.txt
cp ${dir}/${phenoold}.${chr}.${startbp}.${endbp}.uncond.assoc.txt ${pr_dir}/${phenoold}_Neutrophils/1st/${phenoold}.${chr}.${startbp}.${endbp}.uncond.txt
cp ${dir}/${phenoold}.${chr}.${startbp}.${endbp}.uncond.assoc.txt ${pr_dir}/${phenoold}_Monocytes/1st/${phenoold}.${chr}.${startbp}.${endbp}.uncond.txt
rm ${dir}/${phenoold}.${chr}.${startbp}.${endbp}.uncond.txt_head
rm ${dir}/${pheno}.${chr}.${startbp}.${endbp}.uncond.assoc.${plinktest}.ADD.nNA
rm ${dir}/${pheno}.${chr}.${startbp}.${endbp}.uncond.assoc.${plinktest}
rm ${dir}/${pheno}.${chr}.${startbp}.${endbp}.uncond.assoc.${plinktest}.ADD 
rm ${dir}/${pheno}.${chr}.${startbp}.${endbp}.uncond.log


echo $pheno
echo $chr
echo $startbp
echo $endbp

