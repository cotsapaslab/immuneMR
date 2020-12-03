#! /bin/bash


module load PLINK/1.90-beta4.6
module load R

WORK_DIR=$1
PR_DIR=$2
gen1=$3
gen2=$4



# Get overlapping SNPs, remove ambiguous SNPs and Duplicates, Flip
Rscript ${WORK_DIR}/Scripts/getSNPs.R ${gen1}.bim ${gen2}.bim ${PR_DIR}
plink --bfile ${gen1} --extract ${PR_DIR}/keep1.txt --make-bed --out ${PR_DIR}/gen1
plink --bfile ${PR_DIR}/gen1 --exclude ${PR_DIR}/remove_dup1.txt --make-bed --out ${PR_DIR}/gen1
plink --bfile ${PR_DIR}/gen1 --exclude ${PR_DIR}/remove_amb1.txt --make-bed --out ${PR_DIR}/gen1
plink --bfile ${PR_DIR}/gen1 --exclude ${PR_DIR}/remove_diff1.txt --make-bed --out ${PR_DIR}/gen1
plink --bfile ${PR_DIR}/gen1 --flip ${PR_DIR}/flip1.txt --make-bed --out ${PR_DIR}/gen1
plink --bfile ${PR_DIR}/gen1 --update-name ${PR_DIR}/updatenames.txt --make-bed --out ${PR_DIR}/gen1
plink --bfile ${PR_DIR}/gen1 --update-name ${PR_DIR}/updatenamesflips.txt --make-bed --out ${PR_DIR}/gen1
plink --bfile ${gen2} --extract ${PR_DIR}/keep2.txt --make-bed --out ${PR_DIR}/gen2
plink --bfile ${PR_DIR}/gen2 --exclude ${PR_DIR}/remove_dup2.txt --make-bed --out ${PR_DIR}/gen2
plink --bfile ${PR_DIR}/gen2 --exclude ${PR_DIR}/remove_amb2.txt --make-bed --out ${PR_DIR}/gen2
plink --bfile ${PR_DIR}/gen2 --exclude ${PR_DIR}/remove_diff2.txt --make-bed --out ${PR_DIR}/gen2

# Merge datasets to check for duplicates and relatives
plink --bfile ${PR_DIR}/gen1 --bmerge ${PR_DIR}/gen2 --make-bed --out ${PR_DIR}/merged

# Make sure that there is no sample overlap between the data sets
plink --bfile ${PR_DIR}/merged --geno 0.02 --hwe 1e-3 --maf 0.05 --indep-pairwise 200 100 0.2 --out ${PR_DIR}/merged_prunedREL
plink --bfile ${PR_DIR}/merged --extract ${PR_DIR}/merged_prunedREL.prune.in --genome --min 0.0625 --out ${PR_DIR}/merged_genome
plink --bfile ${PR_DIR}/merged --missing --out ${PR_DIR}/merged
Rscript ${WORK_DIR}/Scripts/dupl_and_relatives.R ${PR_DIR}/merged
##::::::::::::::::::::::::::::::::::
# check manually and eventually remove duplicates or relatives (not necessary here)
##::::::::::::::::::::::::::::::::::

