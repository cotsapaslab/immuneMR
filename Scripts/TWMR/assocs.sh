

module load PLINK/1.90-beta4.6
module load R

SNPlist=$1
geneslist=$2
snpslist=$3
DS=$4
gen1=$5
pthres=$6
ImmPhen=$7
ImmPhenDS=$8
ImmPhencov=$9
JLIM=${10}
geno=${11}
WORK_DIR=${12}
PR_DIR=${13}


windows="${WORK_DIR}/${DS}/Data/Pheno/PhenoList_all_windows.txt"
pheno="${WORK_DIR}/${DS}/Data/Pheno/phenotypes_transformed.txt"
covars="${WORK_DIR}/${DS}/Data/Pheno/covariates_selected.txt"
genoImmPhen="${WORK_DIR}/MIP/Data/Genotypes/checkQC/LabExMI_imputation_816x5699237_noMiss_noHet_noRel_noOutl_HWE_MAF_QCd_TWMR"


readarray snps < $SNPlist
numsnps=${#snps[@]}


# Get all cis Genes for these SNPs 

rm $geneslist
touch $geneslist

for (( i=0; i<$numsnps; i++ ))
do
info="${snps[$i]}"
info=$(echo ${info}|tr -d '\n')

OFS=$IFS
IFS=$' '
set -- $info
array=( $@ )
chr=${array[0]}
startbp=${array[3]}
endbp=${array[4]}
IFS=$OFS

#echo "${chr}---${startbp}---${endbp}"

awk -v chr=$chr -v startbp=$startbp -v endbp=$endbp '{ if (NR > 1 && $3 == chr && (($4 > startbp && $4 < endpb) || ($5 > startbp && $5 < endbp))) print $1"\t"$3"\t"$6"\t"$7}' $windows >> $geneslist

done

sort $geneslist | uniq > ${geneslist}_2
mv ${geneslist}_2 $geneslist


# Get all Genes for which at least on of the SNPs is an eQTL

readarray genes < $geneslist
numgenes=${#genes[@]}

rm ${geneslist}_2
touch ${geneslist}_2

for (( i=0; i<$numgenes; i++ ))
do

geneinfo="${genes[$i]}"
geneinfo=$(echo ${geneinfo}|tr -d '\n')
OFS=$IFS
IFS=$' '
set -- $geneinfo
array=( $@ )
gene=${array[0]}
IFS=$OFS

plink --bfile $geno --extract $SNPlist --linear sex --pheno $pheno --pheno-name $gene --covar $covars --out ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp
awk '{if (NR == 1 || $5 == "ADD") print $0}' ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.assoc.linear > ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.assoc.linear.ADD
cat ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.assoc.linear.ADD
mv ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.assoc.linear.ADD ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.assoc.linear
LC_ALL=C sort -k9 -g ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.assoc.linear >  ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.assoc.linear.sort
lead=`awk -v p=$pthres '{ if(NR == 2 && $9<p) print $2}' ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.assoc.linear.sort`
if [ "$lead" != "" ]
then
	if [ "$gene" != "$gen1" ]
	then
	printf "${geneinfo}\n" >> ${geneslist}_2
	fi
fi
rm ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.*
done

mv ${geneslist}_2 ${geneslist}

# Find all other independent SNPs associated  with these Genes


readarray genes < $geneslist
numgenes=${#genes[@]}

awk '{ print $1"\t"$2"\t"$3 }' $SNPlist > $snpslist  


for (( i=0; i<$numgenes; i++ ))
do

geneinfo="${genes[$i]}"
geneinfo=$(echo ${geneinfo}|tr -d '\n')
OFS=$IFS
IFS=$' '
set -- $geneinfo
array=( $@ )
gene=${array[0]}
chr=${array[1]}
startbp=${array[2]}
endbp=${array[3]}
IFS=$OFS

echo "${gene}---${chr}---${startbp}---${endbp}---"

plink --bfile $geno --extract ${genoImmPhen}.snps --linear sex --chr $chr --from-bp $startbp --to-bp $endbp --pheno $pheno --pheno-name $gene --covar $covars --out ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp
awk '{if (NR == 1 || $5 == "ADD") print $0}' ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.assoc.linear > ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.assoc.linear.ADD
bash ${WORK_DIR}/Scripts/TWMR/conditional.sh ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.assoc.linear.ADD $pthres $geno $pheno $gene $covars ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.assoc.linear.indSNPs
echo "--------------"
cat ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.assoc.linear.indSNPs
echo "---------------"
echo ""

cat $snpslist ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.assoc.linear.indSNPs > ${snpslist}_2
mv ${snpslist}_2 $snpslist
rm ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.* 

done

sort $snpslist | uniq > ${snpslist}_2
mv ${snpslist}_2 $snpslist

# Pruning of the final list of genes and create LD file

plink --bfile $geno --extract $snpslist --make-bed --out ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_tmp
plink --bfile ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_tmp --r2 --matrix --out ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_tmp
Rscript ${WORK_DIR}/Scripts/TWMR/getindSNPs.R ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_tmp.ld ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_tmp.bim $SNPlist $snpslist
plink --bfile $geno --extract $snpslist --make-bed --out ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_tmp
plink --bfile ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_tmp --r2 --matrix --out ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_tmp
mv ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_tmp.ld ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.ld
rm ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_tmp.*


# Get assocs of these SNPs to all genes

plink --bfile $geno --extract $snpslist --linear sex --pheno $pheno --pheno-name $gen1 --covar $covars --ci 0.95 --out ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gen1}_tmp
awk '{ if (NR == 1 || $5 == "ADD") print $0}' ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gen1}_tmp.assoc.linear > ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gen1}_tmp.assoc.linear.ADD
mv ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gen1}_tmp.assoc.linear.ADD ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gen1}_tmp.assoc.linear
printf "GENES\t${gen1}\n" > ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gen1}_tmp.header
awk '{if (NR > 1) print $2"\t"$7}' ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gen1}_tmp.assoc.linear > ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gen1}_tmp.restB
awk '{if (NR > 1) print $2"\t"$8}' ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gen1}_tmp.assoc.linear > ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gen1}_tmp.restS
awk '{if (NR > 1) print $2"\t"$12}' ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gen1}_tmp.assoc.linear > ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gen1}_tmp.restP
cat ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gen1}_tmp.header ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gen1}_tmp.restB > ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.matrix
cat ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gen1}_tmp.header ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gen1}_tmp.restP > ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.P
cat ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gen1}_tmp.header ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gen1}_tmp.restS > ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.SE

rm ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gen1}_tmp*


readarray genes < $geneslist
numgenes=${#genes[@]}

for (( i=0; i<$numgenes; i++ ))
do

geneinfo="${genes[$i]}"
geneinfo=$(echo ${geneinfo}|tr -d '\n')
OFS=$IFS
IFS=$' '
set -- $geneinfo
array=( $@ )
gene=${array[0]}
IFS=$OFS

plink --bfile $geno --extract $snpslist --linear sex --pheno $pheno --pheno-name $gene --covar $covars --ci 0.95 --out ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp
awk '{ if (NR == 1 || $5 == "ADD") print $0}' ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.assoc.linear > ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.assoc.linear.ADD
mv ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.assoc.linear.ADD ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.assoc.linear
LC_ALL=C sort -k12 -g ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.assoc.linear > ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.assoc.linear.sort
lead=`awk -v p=$pthres '{ if(NR == 2 && $12<p) print $2}' ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.assoc.linear.sort`
if [ "$lead" != "" ]
then
printf "${gene}\n" > ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.header
awk '{if (NR > 1) print $7}' ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.assoc.linear > ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.restB
awk '{if (NR > 1) print $8}' ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.assoc.linear > ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.restS
awk '{if (NR > 1) print $12}' ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.assoc.linear > ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.restP
cat ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.header ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.restB > ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.addB
cat ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.header ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.restS > ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.addS
cat ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.header ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.restP > ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.addP
paste -d"\t" ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.matrix ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.addB > ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.matrix_2
mv ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.matrix_2 ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.matrix
paste -d"\t" ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.P ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.addP > ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.P_2
mv ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.P_2 ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.P
paste -d"\t" ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.SE ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp.addS > ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.SE_2
mv ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.SE_2 ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.SE


fi

rm ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${gene}_tmp* 

done

# Calculate BETAs for immune phenotye

if [[ $ImmPhen == *"binary"* ]]
then
plink --bfile $genoImmPhen --extract $snpslist --logistic sex --pheno $ImmPhenDS --pheno-name $ImmPhen --covar $ImmPhencov --out ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_tmp
Rscript ${WORK_DIR}/Scripts/TWMR/ORtoBETA.R ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_tmp.assoc.logistic ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_tmp.assoc.linear
else
plink --bfile $genoImmPhen --extract $snpslist --linear sex --pheno $ImmPhenDS --pheno-name $ImmPhen --covar $ImmPhencov --out ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_tmp
fi
awk '{ if (NR == 1 || $5 == "ADD") print $0}' ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_tmp.assoc.linear > ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_tmp.assoc.linear.ADD
mv ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_tmp.assoc.linear.ADD ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_tmp.assoc.linear
printf "BETA_GWAS\n" > ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_tmp.header
awk '{if (NR > 1) print $7}' ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_tmp.assoc.linear > ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_tmp.rest
cat ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_tmp.header ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_tmp.rest > ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_tmp.add
paste -d"\t" ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.matrix ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_tmp.add > ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.matrix_2
mv ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.matrix_2 ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.matrix
rm ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_tmp*


# Remove correlated genes
cp ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.matrix ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${JLIM}.allgenes
Rscript ${WORK_DIR}/Scripts/TWMR/remcorrGenes_byexpr.R ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.matrix $pheno

mv ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.matrix ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${JLIM}.matrix
mv ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.ld ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${JLIM}.ld
mv ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.P ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${JLIM}.P
mv ${PR_DIR}/${DS}_${gen1}_${ImmPhen}.SE ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${JLIM}.SE


rm ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${JLIM}.SE
rm ${PR_DIR}/${DS}_${gen1}_${ImmPhen}_${JLIM}.P

