#! /bin/bash

#SBATCH --job-name=GWAS_getleadSNPs
#SBATCH --output=/home/cg859/scratch60/Logs/GWAS_leads_%A_%a_log.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --array=1-164
#SBATCH --partition=general


module load PLINK/1.90-beta4.6

plinkfile=$1
folder=$2
phenolist=$3

readarray traits < $phenolist

index=$((SLURM_ARRAY_TASK_ID-1))
trait="${traits[$index]}"
trait=$(echo $trait|tr -d '\n')

if [[ $trait == *"binary"* ]]
then
plinktest="logistic"
else
plinktest="linear"
fi

assocfile=${folder}/${trait}/${trait}.${plinktest}

rm -f ${assocfile}_leads
touch ${assocfile}_leads
printf "CHR SNP BP P\n" >> ${assocfile}_leads
i=1

lead="Lead"
if [-f $assocfile ];then
gunzip ${assocfile}.gz
awk '{ if (NR == 1 || $5 == "ADD") print $0}' $assocfile > ${assocfile}_ADD
else
cp ${assocfile}.ADD ${assocfile}_ADD
fi

while [ "$lead" != "" ]
do
plink --bfile $plinkfile --clump ${assocfile}_ADD --clump-p1 0.00001 --clump-r2 0.2 --clump-kb 500 --out ${assocfile}_clump
j=$((i+1))
head -n $j ${assocfile}_clump.clumped > ${assocfile}.clumped
readarray lines < ${assocfile}.clumped
line="${lines[$i]}"
line=$(echo $line|tr -d '\n')
OFS=$IFS
IFS=$' '
set -- $line
array=( $@ )
lead=${array[2]}
chr=${array[0]}
bp=${array[3]}
pval=${array[4]}
IFS=$OFS
startbp=$((bp-100000))
endbp=$((bp+100000))
awk -v chr="$chr" -v sbp="$startbp" -v ebp="$endbp" -v lead="$lead" '{ if (NR == 1 || $2 == lead || $1 != chr || $3 < sbp || $3 > ebp) print $0}' ${assocfile}_ADD > ${assocfile}_ADD_2
mv ${assocfile}_ADD_2 ${assocfile}_ADD
i=$((i+1))
printf "${chr} ${lead} ${bp} ${pval}\n" >> ${assocfile}_leads
rm ${assocfile}_clump.clumped 
rm ${assocfile}.clumped
rm ${assocfile}_clump.log
done
awk '{ if ($2 != "" && $2 != "Lead" && $1 != "") print $0}' ${assocfile}_leads > ${folder}/${trait}/${trait}.leads
rm ${assocfile}_leads
rm ${assocfile}_ADD





