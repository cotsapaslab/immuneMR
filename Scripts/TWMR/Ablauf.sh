#! /bin/bash

#SBATCH --job-name=TWMR
#SBATCH --output=/home/cg859/scratch60/Logs/TWMR_%A_%a_log.txt
#SBATCH --partition=general
#SBATCH --cpus-per-task=1
#SBATCH --time=7:00:00
#SBATCH --mem-per-cpu=5000
#SBATCH --array=2-16653


module load PLINK/1.90-beta4.6
module load R

liste=$1
pthres=$2
WORK_DIR=$3
PR_DIR=$4
genfileImmPhen=$5

i=$((SLURM_ARRAY_TASK_ID))
readarray MRlist < $liste
#numlist=${#MRlist[@]} 

#for (( i=1; i<$numlist; i++ )) 
#do
 
info="${MRlist[$i]}"
info=$(echo $info|tr -d '\n')

OFS=$IFS
IFS=$' '
set -- $info
array=( $@ )
DSname=${array[1]}
gen1=${array[2]}
ImmPhen=${array[0]}
if [[ $ImmPhen != *"binary"* ]]
then
ImmPhen="${ImmPhen}_norm"
fi
JLIM=${array[12]}
IFS=$OFS

	if [ "$DSname" == "Tcells" ]
        then
	DS="BP"
        fi

        if [ "$DSname" == "Monocytes" ]
        then
	DS="BP_Mono"
        fi

        if [ "$DSname" == "Neutrophils" ]
        then
	DS="BP_Neutro"
        fi
 

#echo "${gen1}----${DSname}----${ImmPhen}----${DS}"

assocfile="${WORK_DIR}/${DS}/cis_eQTL/cis${gen1}.assoc.linear.ADD.sort"
genfile="${WORK_DIR}/BP/Data/Geno/checkQC/all_noMiss_noHet_noRel_noOutl_TWMR"
ImmPhenDS="${WORK_DIR}/MIP/Data/Pheno/ImmPhen/phenotypes_transformed.txt"
ImmPhencov="${WORK_DIR}/MIP/Data/Pheno/ImmPhen/covariates_selected.txt"
genfileImmPhen="${WORK_DIR}/MIP/Data/Genotypes/checkQC/LabExMI_imputation_816x5699237_noMiss_noHet_noRel_noOutl_HWE_MAF_QCd_TWMR"
pheno="${WORK_DIR}/${DS}/Data/Pheno/phenotypes_transformed.txt"
covars="${WORK_DIR}/${DS}/Data/Pheno/covariates_selected.txt"
windows="${WORK_DIR}/${DS}/Data/Pheno/PhenoList_all_windows.txt"

# Run cis eQTL analysis on Gene
chr=`awk -v g="$gen1" '{ if ($1 == g) print $3}' $windows`
startbp=`awk -v g="$gen1" '{ if ($1 == g) print $6}' $windows`
endbp=`awk -v g="$gen1" '{ if ($1 == g) print $7}' $windows`
plink --bfile $genfile --extract ${genfileImmPhen}.snps --chr $chr --from-bp $startbp --to-bp $endbp --linear sex --pheno $pheno --pheno-name $gen1 --covar $covars --out ${PR_DIR}/${DS}_${gen1}_${ImmPhen}
awk '{ if (NR == 1 || $5 == "ADD") print $0 }' ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.assoc.linear > ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.assoc.linear.ADD

# Get number of associated SNPs to gene of interest
plink --bfile $genfile --extract ${genfileImmPhen}.snps --make-bed --out ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_shared
bash ${WORK_DIR}/Scripts/TWMR/conditional.sh ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.assoc.linear.ADD $pthres ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_shared $pheno $gen1 $covars ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.SNPs
awk '{ print $1"\t"$2"\t"$3"\t"($2-1000000)"\t"($2+1000000) }' ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.SNPs > ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.SNPs_2
mv ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.SNPs_2 ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.SNPs
 
# Create TWMR files 
bash ${WORK_DIR}/Scripts/TWMR/assocs.sh ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.SNPs ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.Genes ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.allSNPs ${DS} ${gen1} $pthres $ImmPhen $ImmPhenDS $ImmPhencov $JLIM ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_shared ${WORK_DIR} ${PR_DIR}
rm ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.*
rm ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_shared*

# Switch BETAs for "mismatch" SNPs
Rscript ${WORK_DIR}/Scripts/TWMR/switchbetas.R ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${JLIM}.matrix ${PR_DIR}/../mismatch_snps.txt

# Run TWMR
cd ${PR_DIR}
Rscript ${WORK_DIR}/Scripts/TWMR/TWMR.R ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${JLIM} 816 197 

#fi

