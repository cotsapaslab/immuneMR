#! /bin/bash


#SBATCH --job-name=PRSice
#SBATCH --output=${WORK_DIR}/Logs/PRS_run_%A_%a_log.txt
#SBATCH --partition=general
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=5000
#SBATCH --array=0-1665

module load PLINK/1.90-beta4.6
module load R

SECONDS=0

startindex=$((SLURM_ARRAY_TASK_ID * 10))
startindex=$((startindex + 1))
endindex=$((startindex + 10))

jlimfile=$1
readarray lines < ${jlimfile}

WORK_DIR=$2
PR_DIR=$3

for (( cpi=$startindex; cpi<$endindex; cpi++ ))
do

index=$cpi
#index=$SLURM_ARRAY_TASK_ID

   line="${lines[$index]}"
   line=$(echo ${line}|tr -d '\n')


          	filename=${line}
                OFS=$IFS
                IFS=$' '
                set -- $filename
                array=( $@ )
                maintr=${array[0]}
		if [[ $maintr != *"binary"* ]]
		then
		maintr="${maintr}_norm"
		binarytrait="F"
		else
		binarytrait="T"
		fi
                sectr=${array[2]}
                dataset=${array[1]}
                idxSNP=${array[4]}
                idxBP=${array[5]}
                chr=${array[3]}
		condlead=${array[9]}
		IFS=$OFS 

origgen="${WORK_DIR}/MIP/Data/Genotypes/checkQC/LabExMI_imputation_816x5699237_noMiss_noHet_noRel_noOutl"

pheno2="${WORK_DIR}/MIP/Data/Pheno/ImmPhen/phenotypes_transformed.txt"
cov2="${WORK_DIR}/MIP/Data/Pheno/ImmPhen/covariates_selected.txt"

dataset1=$dataset
phenotype1=$sectr
phenotype2=$maintr
jlimSNP=$idxSNP
indexbp=$idxBP
chr=$chr


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
	


	locusstart=$((indexbp-100000))
	locusend=$((indexbp+100000))
	locus="${chr}.${locusstart}.${locusend}"

	cp ${PR_DIR}/${dataset1}_${phenotype1}_cond${condlead}.assoc.linear.gz ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}_${locus}.assoc.linear_condJLIM.gz
	gunzip ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}_${locus}.assoc.linear_condJLIM.gz
	plink --bfile ${PR_DIR}/gen2 --make-bed --out ${PR_DIR}/gen2_${dataset1}_${phenotype1}_${phenotype2}_${locus}_condJLIM

Rscript ${WORK_DIR}/Scripts/PRS/PRSice.R --dir . --prsice ${WORK_DIR}/Scripts/PRS/PRSice_linux --base ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}_${locus}.assoc.linear_condJLIM --target ${PR_DIR}/gen2_${dataset1}_${phenotype1}_${phenotype2}_${locus}_condJLIM --pheno ${pheno2} --pheno-col ${phenotype2} --cov ${cov2} --A1 A1 --A2 A2 --model add --score sum --thread 1 --stat BETA --beta --binary-target ${binarytrait} --bar-levels 0.01,0.001,0.0001,0.00005,0.00001,0.000005,0.000001,0.0000005,0.0000001,0.00000005 --fastscore --no-full --out ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}_${locus}_2pheno_condJLIMlocus
awk '{ if (NR != 1) print $2 }' ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}_${locus}_2pheno_condJLIMlocus.prsice > ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}_${locus}_2pheno_condJLIMlocus.thresh

readarray pthreshs < ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}_${locus}_2pheno_condJLIMlocus.thresh

for pthres in "5e-08" "1e-07" "5e-07" "1e-06" "5e-06" "1e-05" "5e-05" "0.0001" "0.001" "0.01"
do

found="0"
for i in "${pthreshs[@]}"
do
        akt=$i
        akt=$(echo $akt|tr -d '\n')
    if [ "$akt" == "$pthres" ] ; then
        found="1"
        echo $akt
    fi
done

if [ "$found" == "0" ] ; then

# Calculate R2 for second phenotype
Rscript ${WORK_DIR}/Scripts/PRS/PRSice.R --dir . --prsice ${WORK_DIR}/Scripts/PRS/PRSice_linux --base ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}_${locus}.assoc.linear_condJLIM --target ${PR_DIR}/gen2_${dataset1}_${phenotype1}_${phenotype2}_${locus}_condJLIM --pheno ${pheno2} --pheno-col ${phenotype2} --cov ${cov2} --A1 A1 --A2 A2 --model add --score sum --thread 1 --stat BETA --beta --binary-target ${binarytrait} --bar-levels ${pthres} --fastscore --no-full --lower ${pthres} --upper ${pthres} --print-snp --print-snp --out ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}_${locus}_2pheno_condJLIMlocus_p_${pthres}
# --print-snp


if [ -f ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}_${locus}_2pheno_condJLIMlocus_p_${pthres}.prsice ]
then

### what do I need this file for??
#mv ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}_${locus}_2pheno_condJLIMlocus_p_${pthres}.snp ${PR_DIR}/snplist_best_${dataset1}_${phenotype1}_${phenotype2}_${locus}_2pheno_condJLIMlocus_p_${pthres}.txt

awk '{if (NR != 1) print $0}' ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}_${locus}_2pheno_condJLIMlocus_p_${pthres}.prsice > ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}_${locus}_2pheno_condJLIMlocus_p_${pthres}.prsice_2
cat ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}_${locus}_2pheno_condJLIMlocus.prsice ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}_${locus}_2pheno_condJLIMlocus_p_${pthres}.prsice_2 > ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}_${locus}_2pheno_condJLIMlocus.prsice_2
awk '{ if (NR == 1 || $1 != "Set") print $0 }' ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}_${locus}_2pheno_condJLIMlocus.prsice_2 > ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}_${locus}_2pheno_condJLIMlocus.prsice
rm ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}_${locus}_2pheno_condJLIMlocus.prsice_2
rm ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}_${locus}_2pheno_condJLIMlocus_p_${pthres}.prsice_2

fi

rm -f ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}_${locus}_2pheno_condJLIMlocus_p*.prsice
rm -f ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}_${locus}_2pheno_condJLIMlocus_p*.png
rm -f ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}_${locus}_2pheno_condJLIMlocus_p*.best
rm -f ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}_${locus}_2pheno_condJLIMlocus_p*.log
rm -f ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}_${locus}_2pheno_condJLIMlocus_p*.summary

fi

done

LC_ALL=C sort -k2 -g ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}_${locus}_2pheno_condJLIMlocus.prsice > ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}_${locus}_2pheno_condJLIMlocus.prsice_2
mv ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}_${locus}_2pheno_condJLIMlocus.prsice_2 ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}_${locus}_2pheno_condJLIMlocus.prsice

rm ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}_${locus}.assoc.linear_condJLIM
rm ${PR_DIR}/gen2_${dataset1}_${phenotype1}_${phenotype2}_${locus}_condJLIM.log
rm ${PR_DIR}/gen2_${dataset1}_${phenotype1}_${phenotype2}_${locus}_condJLIM.bim
rm ${PR_DIR}/gen2_${dataset1}_${phenotype1}_${phenotype2}_${locus}_condJLIM.bed
rm ${PR_DIR}/gen2_${dataset1}_${phenotype1}_${phenotype2}_${locus}_condJLIM.fam
rm ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}.${locus}.JLIMlocus.fam
rm ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}.${locus}.JLIMlocus.log
rm ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}.${locus}.JLIMlocus.bim
rm ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}.${locus}.JLIMlocus.bed
rm ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}.${locus}.JLIMlocus.*
rm ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}_${locus}_2pheno_condJLIMlocus.summary
rm ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}_${locus}_2pheno_condJLIMlocus*.png
rm ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}_${locus}_2pheno_condJLIMlocus.best
rm ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}_${locus}_2pheno_condJLIMlocus.log
rm ${PR_DIR}/${dataset1}_${phenotype1}_${phenotype2}_${locus}_2pheno_condJLIMlocus.thresh

done

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."


