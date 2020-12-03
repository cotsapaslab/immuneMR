#!/bin/bash

#SBATCH --partition=general
#SBATCH --time=1:00:00
#SBATCH --job-name=indep_1st
#SBATCH --output=/home/cg859/scratch60/Logs/PRS_indep2nd_%A_%a_log
#SBATCH --array=1-6

#12928 Tcells
#9990 Neutrophils
#12201 Monocytes

plinkfile=$2
phenofile=$3
covarfile=$4
WD=$5
regdir=$6
readarray lines < $1

module load PLINK/1.90-beta5.3

index=$((SLURM_ARRAY_TASK_ID-1))

line="${lines[$index]}"
line=$(echo $line|tr -d '\n')
echo $line
OFS=$IFS
IFS=$' '
set -- $line
array=( $@ )
IFS=$OFS

gene=${array[0]}
chr=${array[1]}
bp=${array[2]}
startbp=$((bp-1200000))
endbp=$((bp+1200000))

echo $chr
echo $startbp
echo $endbp
echo $gene

if (( startbp < 0 ))
then
startbp=0
fi

bash ${WD}/Scripts/conditional_2nd.sh $plinkfile $phenofile $covarfile $gene $chr $startbp $endbp $regdir $WD


# Run unconditional analysis
plinktest="linear"
plink --bfile ${plinkfile} --${plinktest} sex --pheno ${phenofile} --pheno-name ${gene} --covar ${covarfile} --ci 0.95 --chr $chr --from-bp $startbp --to-bp $endbp --out ${regdir}/${gene}.uncond
grep ADD ${regdir}/${gene}.uncond.assoc.${plinktest} > ${regdir}/${gene}.uncond.assoc.${plinktest}.ADD
grep -v NA ${regdir}/${gene}.uncond.assoc.${plinktest}.ADD > ${regdir}/${gene}.uncond.assoc.${plinktest}.ADD.nNA
printf "CHR\tSNP\tBP\tA1\tTEST\tNMISS\tBETA\tSE\tLxy\tHxy\tSTAT\tP\n" > ${regdir}/${gene}.uncond.txt_head
cat ${regdir}/${gene}.uncond.txt_head ${regdir}/${gene}.uncond.assoc.${plinktest}.ADD.nNA > ${regdir}/${gene}.uncond.assoc.txt
rm ${regdir}/${gene}.uncond.txt_head
rm ${regdir}/${gene}.uncond.assoc.${plinktest}.ADD.nNA
rm ${regdir}/${gene}.uncond.assoc.${plinktest}
rm ${regdir}/${gene}.uncond.assoc.${plinktest}.ADD
rm ${regdir}/${gene}.uncond.log 

# Run n-1 conditional analyses
numind=$(wc -l < "${regdir}/${gene}.conds.txt")

if [ "$numind" != "1" ]; then
for (( i=1; i<=$numind; i++ ))
do
awk -v i="$i" '{ if (NR != i) print $0}' ${regdir}/${gene}.conds.txt > ${regdir}/${gene}.conds.txt_${i}
plink --bfile ${plinkfile} --${plinktest} sex --pheno ${phenofile} --pheno-name ${gene} --covar ${covarfile} --ci 0.95 --chr $chr --from-bp $startbp --to-bp $endbp --condition-list ${regdir}/${gene}.conds.txt_${i} --out ${regdir}/${gene}.cond${i}
grep ADD ${regdir}/${gene}.cond${i}.assoc.${plinktest} > ${regdir}/${gene}.cond${i}.assoc.${plinktest}.ADD
grep -v NA ${regdir}/${gene}.cond${i}.assoc.${plinktest}.ADD > ${regdir}/${gene}.cond${i}.assoc.${plinktest}.ADD.nNA
rm ${regdir}/${gene}.cond${i}.assoc.${plinktest}.ADD
rm ${regdir}/${gene}.cond${i}.assoc.${plinktest}
rm ${regdir}/${gene}.cond${i}.log

printf "CHR\tSNP\tBP\tA1\tTEST\tNMISS\tBETA\tSE\tLxy\tHxy\tSTAT\tP\n" > ${regdir}/${gene}.cond${i}.txt_head
cat ${regdir}/${gene}.cond${i}.txt_head ${regdir}/${gene}.cond${i}.assoc.${plinktest}.ADD.nNA >  ${regdir}/${gene}.cond${i}.assoc.txt
rm ${regdir}/${gene}.cond${i}.assoc.${plinktest}.ADD.nNA
rm ${regdir}/${gene}.cond${i}.txt_head

done
fi

echo "DONE"
