#! /bin/bash

module load PLINK/1.90-beta4.6
module load R



ImmPhen=$1
Gene=$2
DS=$3
cond1st=$4
cond2nd=$5
chr=$6
startbp=$7
endbp=$8
firsttraitp="$9"
firsttraitn="${10}"
phenoG="${11}"
PR_DIR="${12}"
WORK_DIR="${13}"
GWAS_DIR="${14}"
gen1file="${15}"

ass1file="${GWAS_DIR}/${ImmPhen}.${chr}.${startbp}.${endbp}.${cond1st}.assoc.txt"
ass2file="${PR_DIR}/${DS}_assoc/${Gene}.${cond2nd}.${chr}.${startbp}.${endbp}.assoc.txt"
pheno=$Gene
ldfile="${WORK_DIR}/PRS-Project/ld0/locus.${chr}.${startbp}.${endbp}.txt.gz"
gunzip $ldfile
ldfile="${WORK_DIR}/PRS-Project/ld0/locus.${chr}.${startbp}.${endbp}.txt"


# Add A2 to ass1 file
Rscript ${WORK_DIR}/Scripts/addA2_plots.R ${WORK_DIR}/MIP/Data/Genotypes/checkQC/LabExMI_imputation_816x5699237_noMiss_noHet_noRel_noOutl_HWE_MAF_QCd.bim $ass1file ${ass1file}_A2
ass1file="${ass1file}_A2"

# Create common lead variant file
Rscript ${WORK_DIR}/Scripts/commonleadvariant.R $ass1file $ass2file ${PR_DIR}/Figure1_lead.txt
leadfile="${PR_DIR}/Figure1_lead.txt"
leadsnp=$(head -n 1 ${leadfile})

# Create LD file (in phase information) for lead SNP
if [ ! -f ${PR_DIR}/${ImmPhen}_${DS}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${Gene}.${cond2nd}/2nd/locus.${chr}.${startbp}.${endbp}/2nd.ped.gz ]; then
	if [ ! -f ${PR_DIR}/${ImmPhen}_${DS}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${Gene}.${cond2nd}/2nd/locus.${chr}.${startbp}.${endbp}/2nd.ped ]; then
	cp ${PR_DIR}/${ImmPhen}_${DS}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${Gene}.${cond2nd}/2nd/locus.${chr}.${startbp}.${endbp}/*.ped.gz ${PR_DIR}/${ImmPhen}_${DS}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${Gene}.${cond2nd}/2nd/locus.${chr}.${startbp}.${endbp}/2nd.ped.gz
	gunzip ${PR_DIR}/${ImmPhen}_${DS}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${Gene}.${cond2nd}/2nd/locus.${chr}.${startbp}.${endbp}/2nd.ped.gz
	fi
fi

if [ ! -f ${PR_DIR}/${ImmPhen}_${DS}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${Gene}.${cond2nd}/2nd/locus.${chr}.${startbp}.${endbp}/2nd.ped ]; then
gunzip ${PR_DIR}/${ImmPhen}_${DS}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${Gene}.${cond2nd}/2nd/locus.${chr}.${startbp}.${endbp}/2nd.ped.gz
fi

plink --file ${PR_DIR}/${ImmPhen}_${DS}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${Gene}.${cond2nd}/2nd/locus.${chr}.${startbp}.${endbp}/2nd --r in-phase --ld-snp ${leadsnp} --inter-chr --out ${PR_DIR}/Figure1_phase
ldplinkfile=${PR_DIR}/Figure1_phase.ld

# Create 1st trait ped file
plink --bfile $gen1file --chr $chr --from-bp $startbp --to-bp $endbp --recode tab --out ${WORK_DIR}/Figure1_geno1

# Create plot
pedfileG="${PR_DIR}/${ImmPhen}_${DS}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${Gene}.${cond2nd}/2nd/locus.${chr}.${startbp}.${endbp}/2nd.ped"
mapfileG="${PR_DIR}/${ImmPhen}_${DS}/locus.${chr}.${startbp}.${endbp}.${cond1st}/2nd/locus.${chr}.${startbp}.${endbp}.${Gene}.${cond2nd}/2nd/locus.${chr}.${startbp}.${endbp}/2nd.map"
phenofileG=$phenoG
pedfileI="${WORK_DIR}/Figure1_geno1.ped"
mapfileI="${WORK_DIR}/Figure1_geno1.map"
phenofileI="${WORK_DIR}/MIP/Data/Pheno/ImmPhen/phenotypes_transformed.txt"

outputfolder="${WORK_DIR}/Figure_JLIM/${ImmPhen}_${Gene}_${DS}"
mkdir ${WORK_DIR}/Figure_JLIM
mkdir $outputfolder
outputfolder="${outputfolder}/"
Rscript ${WORK_DIR}/Scripts/Submission/plot_Figure_JLIM.R $ass1file $ass2file $Gene $ldfile $ldplinkfile $leadfile $ImmPhen $outputfolder $pedfileG $mapfileG $phenofileG $pedfileI $mapfileI $phenofileI $firsttraitn "$firsttraitp" $DS

gzip ${ldfile}
#rm ${PR_DIR}/Figure1_phase.*
#rm $leadfile






