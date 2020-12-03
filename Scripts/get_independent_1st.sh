#!/bin/bash

#SBATCH --partition=general
#SBATCH --time=72:00:00
#SBATCH --job-name=indep_1st
#SBATCH --output=/home/cg859/scratch60/Logs/PRS_indep1st_%A_%a_log
#SBATCH --array=2-494



DIR=$2
plinkfile=$3
phenofile=$4
covarfile=$5
WD=$6
regdir=$7
readarray lines < $1



index=$((SLURM_ARRAY_TASK_ID-1))
#index=1

line="${lines[$index]}"
line=$(echo $line|tr -d '\n')
echo $line
OFS=$IFS
IFS=$' '
set -- $line
array=( $@ )
IFS=$OFS

firstname=${array[0]}
secondname=${array[1]}

echo $secondname

if [[ $secondname == "Tcells" ]]
then

readarray loci < ${DIR}/${firstname}_Tcells/1st/${firstname}-index.tsv
numloci=${#loci[@]}

if [[ $firstname != *"binary"* ]]
then
firstname=${firstname}_norm
fi


for (( i=1; i<$numloci; i++ ))
do

index2=$i

locus="${loci[$index2]}"
locus=$(echo $locus|tr -d '\n')
OFS=$IFS
IFS=$' '
set -- $locus
array=( $@ )
chr=${array[0]}
startbp=${array[3]}
endbp=${array[4]}
IFS=$OFS

echo $chr
echo $startbp
echo $firstname


bash ${WD}/Scripts/conditional_1st.sh $plinkfile $phenofile $covarfile $firstname $chr $startbp $endbp $regdir $WD


done

fi

#print "DONE"
