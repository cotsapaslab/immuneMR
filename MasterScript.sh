
#####################################################################
# This script lists and describes all the steps perforemed in the project. It is not intended as a stand-alone, unsupervised piece of code - rather, it steps through the various procedures and the single steps should be called one after the other by the user. User input is required at variuos stages. 
#####################################################################


module load PLINK/1.90-beta4.6
module load R
module load Pandoc
module load VCFtools
module load tabix


## Paths to pipeline directories
WORK_DIR="/home/project" # Home directory for the project
# The scripts to run the pipeline need to be in $WORK_DIR/Scripts
PR_DIR="${WORK_DIR}/MIP_BP" # Many of the scripts point to this directory - all results are saved here. 
ImmPhen_DataDir="${WORK_DIR}/MIP/Data/Pheno/ImmPhen" # Directory for the Milieu Intérieur project phenotype data
# LabExMI_covariates.txt needs to be in ${WORK_DIR}/MIP/Data/Pheno/
ImmPhen_DIR="${WORK_DIR}/MIP" # Directory for all data (raw and produced) for the Milieu Intérieur project data set, including GWAS results
DATA_DIR="${WORK_DIR}/BP/Data" # Directory for the BLUEPRINT data set 


## Original data provided by BLUEPRINT and the Milieu Intérieur project needs to be put in the following locations:
# Milieu Intérieur project data:
	
	## Genotype data:
	# LabExMI_imputation_816x5699237.bed in ${WORK_DIR}/MIP/Data/Genotypes/ 
	# LabExMI_imputation_816x5699237.bim in ${WORK_DIR}/MIP/Data/Genotypes/
	# LabExMI_imputation_816x5699237.fam in ${WORK_DIR}/MIP/Data/Genotypes/
	# LabExMI_imputation_816x5699237.tgz in ${WORK_DIR}/MIP/Data/Genotypes/
	
	## Phenotype data:
	# LabExMI_covariates.txt in ${WORK_DIR}/MIP/Data/Pheno/
	# LabExMI_rawfacs.txt in ${WORK_DIR}/MIP/Data/Pheno/ImmPhen/

# BLUEPRINT data: 
	
	## Genotype data:
	#  EGAZ00001235598_release_vcf_06092016_All_chr.BPWP10_13_12_15.vcf in ${WORK_DIR}/BP/Data/Geno

	## Phenotype data: 
	## This data is produced as described in Scripts/GenExp_phenodata.sh and the Gtex-Pipeline, the following results files need to be provided:
	# For T cells: all.expression.bed and all.PEER_covariates.txt in ${WORK_DIR}/BP/Data/Pheno/Gtex/
	#	       phenotypes_transformed.txt in ${WORK_DIR}/BP/Data/Pheno/
	# For monocytes: all.expression.bed and all.PEER_covariates.txt in ${WORK_DIR}/BP_Mono/Data/Pheno/Gtex/
        #              phenotypes_transformed.txt in ${WORK_DIR}/BP_Mono/Data/Pheno/
	# For neutrophils: all.expression.bed and all.PEER_covariates.txt in ${WORK_DIR}/BP_Neutro/Data/Pheno/Gtex/
        #              phenotypes_transformed.txt in ${WORK_DIR}/BP_Neutro/Data/Pheno/	       
	# Metadata of the BLUEPRINT samples can be obtained from the BLUEPRINT consortium (directory called EGAD00001002663 and needs to be in ${WORK_DIR}/BP/Data/


## Software/executables that need to be downloaded and sotred: 
	# PRSice-2: PRSice.R and PRSice_linux in ${WORK_DIR}/Scripts/PRS/
	# JLIM 2.0: download (directory called jlim-master) and put into ${WORK_DIR}/Scripts/JLIM/


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 1. Prepare MIP phenotypes
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

### Phenotype data has been provided by the Milieu Interieur project, the file is called LabExMI_rawfacs.txt and stored in $ImmPhen_DataDir (see above)

# Perform invers-rank transformation, remove two traits (after visual inspection of the phenotype distribution plots) and binarize two traits
Rscript ${WORK_DIR}/Scripts/preparePhenosImmMIP.R ${ImmPhen_DIR}/Data ${WORK_DIR}

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 2. Prepare MIP Genotypes
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   
### Imputed genotype data has been provided by the Milieu Interieur project (plink format), the files are called LabExMI_imputation_816x5699237.bim/.bed./fam and are stored in ${WORK_DIR}/MIP/Data/Genotypes/
Rscript -e "rmarkdown::render('${WORK_DIR}/Scripts/MIP_SampleQC.Rmd', clean=TRUE)" #Datafolder information needs to be changed in the script
plink --bfile ${ImmPhen_DIR}/Data/Genotypes/checkQC/LabExMI_imputation_816x5699237_noMiss_noHet_noRel_noOutl --hwe 0.00001 --maf 0.05 --make-bed --out ${ImmPhen_DIR}/Data/Genotypes/checkQC/LabExMI_imputation_816x5699237_noMiss_noHet_noRel_noOutl_HWE_MAF
plink --bfile ${ImmPhen_DIR}/Data/Genotypes/checkQC/LabExMI_imputation_816x5699237_noMiss_noHet_noRel_noOutl_HWE_MAF --freq --out ${ImmPhen_DIR}/Data/Genotypes/checkQC/LabExMI_imputation_816x5699237_noMiss_noHet_noRel_noOutl_HWE_MAF
Rscript ${WORK_DIR}/Scripts/get_indels.R ${ImmPhen_DIR}/Data/Genotypes/checkQC/LabExMI_imputation_816x5699237_noMiss_noHet_noRel_noOutl_HWE_MAF.bim ${ImmPhen_DIR}/Data/Genotypes/checkQC/indels.txt
Rscript ${WORK_DIR}/Scripts/find_multiallelic.R ${ImmPhen_DIR}/Data/Genotypes/checkQC/LabExMI_imputation_816x5699237_noMiss_noHet_noRel_noOutl_HWE_MAF.bim ${ImmPhen_DIR}/Data/Genotypes/checkQC/multiallelics.txt #no multiallelic SNPs found
plink --bfile ${ImmPhen_DIR}/Data/Genotypes/checkQC/LabExMI_imputation_816x5699237_noMiss_noHet_noRel_noOutl_HWE_MAF --exclude ${ImmPhen_DIR}/Data/Genotypes/checkQC/indels.txt --make-bed --out ${ImmPhen_DIR}/Data/Genotypes/checkQC/LabExMI_imputation_816x5699237_noMiss_noHet_noRel_noOutl_HWE_MAF_QCd

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 3. Prepare MIP Covariates
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# calculate MDS and PCS
bash ${WORK_DIR}/Scripts/MDS_plink.sh ${WORK_DIR}/Scripts/remove_prune.txt ${ImmPhen_DIR}/Data/Genotypes/checkQC/ LabExMI_imputation_816x5699237_noMiss_noHet_noRel_noOutl

# create covariate file (inlucding first 5 PCA, Age, CMV status and tabac use)
Rscript ${WORK_DIR}/Scripts/createCovarMIP.R ${ImmPhen_DIR}/Data/Genotypes/checkQC/LabExMI_imputation_816x5699237_noMiss_noHet_noRel_noOutl_cluster.mds ${ImmPhen_DIR}/Data/Pheno/LabExMI_covariates.txt ${ImmPhen_DIR}/Data/Genotypes/checkQC/LabExMI_imputation_816x5699237_noMiss_noHet_noRel_noOutl_cluster.eigenvec pca ${ImmPhen_DIR}/Data/Pheno/Cell_proportions.txt ${ImmPhen_DIR}/Data/Pheno
awk '{ print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10 }' ${ImmPhen_DIR}/Data/Pheno/covariates.txt > ${ImmPhen_DIR}/Data/Pheno/ImmPhen/covariates_selected.txt

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 4. Prepare BLUEPRINT Phenotypes
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# See  ${WORK_DIR}/Scripts/GenExp_phenodata.sh
# We downloaded FASTQ files and used the GTEx pipeline for RNA-seq alignment, quantification and quality control (see and https://www.gtexportal.org/, Analysis Methods for V8 and https://github.com/broadinstitute/gtex-pipeline).
# The Script ${WORK_DIR}/Scripts/GenExp_phenodata.sh is intended to give an overview on how we applied the GTEx pipeline and can be used for orientation, but is not directly applicable to other data sets and needs user input or small changes at various stages depending on the data. This script was used for the BluePrint T cell data.

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 5. Prepare BLUEPRINT Genotypes
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

### BLUEPRINT genotype data was provided as EGAZ00001235598_release_vcf_06092016_All_chr.BPWP10_13_12_15.vcf.gz and is stored in ${WORK_DIR}/BP/Data/Geno

# Create plink files
gunzip ${DATA_DIR}/Geno/EGAZ00001235598_release_vcf_06092016_All_chr.BPWP10_13_12_15.vcf.gz
plink --vcf ${DATA_DIR}/Geno/EGAZ00001235598_release_vcf_06092016_All_chr.BPWP10_13_12_15.vcf --double-id --make-bed --out ${DATA_DIR}/Geno/all

plink --bfile ${DATA_DIR}/Geno/all --maf 0.05 --hwe 0.00001 --make-bed --out ${DATA_DIR}/Geno/all
#remove indels
echo "." > ${DATA_DIR}/Geno/indelname.txt
plink --bfile ${DATA_DIR}/Geno/all --exclude ${DATA_DIR}/Geno/indelname.txt --make-bed --out ${DATA_DIR}/Geno/all
Rscript ${WORK_DIR}/Scripts/get_indels.R ${DATA_DIR}/Geno/all.bim ${DATA_DIR}/Geno/rem_indels.txt
plink --bfile ${DATA_DIR}/Geno/all --exclude ${DATA_DIR}/Geno/rem_indels.txt --make-bed --out ${DATA_DIR}/Geno/all
rm ${DATA_DIR}/Geno/rem_indels.txt
rm ${DATA_DIR}/Geno/*~
# Add sex information to FAM and create Infofile
Rscript ${WORK_DIR}/Scripts/createInfoFileBP.R ${DATA_DIR}/EGAD00001002663/delimited_maps/Run_Sample_meta_info.map ${DATA_DIR}/Geno/all.fam ${DATA_DIR}/Pheno/infoFile.txt

mkdir ${DATA_DIR}/Geno/checkQC
cp ${WORK_DIR}/Scripts/remove_prune.txt ${DATA_DIR}/Geno/
Rscript -e "rmarkdown::render('${WORK_DIR}/Scripts/BP_SampleQC.Rmd', clean=TRUE)" #Datafolder information needs to be changed in the script

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 6. Prepare BLUEPRINT Covariates
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Calculate MDS and PCA components
bash ${WORK_DIR}/Scripts/MDS_plink.sh ${DATA_DIR}/Geno/remove_prune.txt ${DATA_DIR}/Geno/checkQC/ all_noMiss_noHet_noRel_noOutl

#create covariate-File (first five PCA components, Age, and 30 PEER-factors)
Rscript ${WORK_DIR}/Scripts/createCovarBP.R ${DATA_DIR}/Geno/checkQC/all_noMiss_noHet_noRel_noOutl.fam ${DATA_DIR}/Pheno/infoFile.txt ${DATA_DIR}/Geno/checkQC/all_noMiss_noHet_noRel_noOutl_cluster.mds ${DATA_DIR}/Geno/checkQC/all_noMiss_noHet_noRel_noOutl_cluster.eigenvec pca ${DATA_DIR}/Pheno/covariates.txt

# Select covariates
Rscript ${WORK_DIR}/Scripts/covariatesBP.R ${WORK_DIR}/BP/Data/Pheno/covariates.txt ${WORK_DIR}/BP/Data/Pheno/Gtex/all.PEER_covariates.txt ${WORK_DIR}/BP/Data/Pheno/covariates_selected.txt
Rscript ${WORK_DIR}/Scripts/covariatesBP.R ${WORK_DIR}/BP/Data/Pheno/covariates.txt ${WORK_DIR}/BP_Mono/Data/Pheno/Gtex/all.PEER_covariates.txt ${WORK_DIR}/BP_Mono/Data/Pheno/covariates_selected.txt
Rscript ${WORK_DIR}/Scripts/covariatesBP.R ${WORK_DIR}/BP/Data/Pheno/covariates.txt ${WORK_DIR}/BP_Neutro/Data/Pheno/Gtex/all.PEER_covariates.txt ${WORK_DIR}/BP_Neutro/Data/Pheno/covariates_selected.tx

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 7. Run GWAS on all MIP phenotypes
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Run regression analyses 
mkdir ${ImmPhen_DIR}/GWAS
mkdir ${ImmPhen_DIR}/GWAS/regression
#SLURM
sbatch ${WORK_DIR}/Scripts/plink_regression_MIP.sh ${ImmPhen_DIR}/GWAS ${ImmPhen_DIR}/Data/Genotypes/checkQC/LabExMI_imputation_816x5699237_noMiss_noHet_noRel_noOutl_HWE_MAF_QCd ${ImmPhen_DIR}/Data/Pheno/ImmPhen ${ImmPhen_DIR}/Data/Pheno/ImmPhen/PhenoList.txt
#SLURM
sbatch ${WORK_DIR}/Scripts/plink_combine.sh ${WORK_DIR} ${ImmPhen_DIR}/GWAS ${ImmPhen_DIR}/Data/Pheno/ImmPhen/PhenoList.txt 

# Get lead SNPs
#SLURM
sbatch ${WORK_DIR}/Scripts/getleadSNPs_first.sh ${ImmPhen_DIR}/Data/Genotypes/checkQC/LabExMI_imputation_816x5699237_noMiss_noHet_noRel_noOutl_HWE_MAF_QCd ${ImmPhen_DIR}/GWAS/regression ${ImmPhen_DIR}/Data/Pheno/ImmPhen/PhenoList.txt

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 8. Create List of loci with association
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

### Directory information in this script has to be changed manually
bash ${WORK_DIR}/Scripts/createlistMIP.sh ${ImmPhen_DataDir}/PhenoList.txt ${ImmPhen_DIR}/GWAS ${PR_DIR}/list.txt

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 9. Perform conditional analyses on primary traits (MIP immune phenotypes)
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Get Loci information
#SLURM
sbatch ${WORK_DIR}/Scripts/Loci_1st_part1.sh ${PR_DIR}/list.txt ${WORK_DIR} ${PR_DIR} ${WORK_DIR}/Scripts

# Run conditional analysis on loci to determine independent signals
#SLURM
sbatch ${WORK_DIR}/Scripts/get_independent_1st.sh ${PR_DIR}/list.txt ${PR_DIR} ${ImmPhen_DIR}/Data/Genotypes/checkQC/LabExMI_imputation_816x5699237_noMiss_noHet_noRel_noOutl_HWE_MAF_QCd ${ImmPhen_DataDir}/phenotypes_transformed.txt ${ImmPhen_DataDir}/covariates_selected.txt ${WORK_DIR} ${ImmPhen_DIR}/GWAS/regression 

# Get loci-phenotype combinations with more than one independent signal
wc -l ${ImmPhen_DIR}/GWAS/regression/*conds* > ${PR_DIR}/list_independent_1st.txt
awk '{ if ($1 > 1 && $2 != "total" && $2 != "insgesamt") print $0 }' ${PR_DIR}/list_independent_1st.txt > ${PR_DIR}/list_independent_1st.txt_2
mv ${PR_DIR}/list_independent_1st.txt_2 ${PR_DIR}/list_independent_1st.txt
wc -l ${ImmPhen_DIR}/GWAS/regression/*conds* > ${PR_DIR}/list_single_1st.txt
awk '{ if ($1 == 1 && $2 != "total" && $2 != "insgesamt") print $0 }' ${PR_DIR}/list_single_1st.txt > ${PR_DIR}/list_single_1st.txt_2
mv ${PR_DIR}/list_single_1st.txt_2 ${PR_DIR}/list_single_1st.txt

# Create association stat files conditioning on all combinations of lead SNPs 
#SLURM
sbatch ${WORK_DIR}/Scripts/conditional_assocs_1st.sh ${PR_DIR}/list_independent_1st.txt ${ImmPhen_DIR}/Data/Genotypes/checkQC/LabExMI_imputation_816x5699237_noMiss_noHet_noRel_noOutl_HWE_MAF_QCd ${ImmPhen_DataDir}/phenotypes_transformed.txt ${ImmPhen_DataDir}/covariates_selected.txt ${ImmPhen_DIR}/GWAS/regression ${PR_DIR} 

# Create unconditional association stat files for loci with a single signal
#SLURM
sbatch ${WORK_DIR}/Scripts/unconditional_assocs_1st.sh ${PR_DIR}/list_single_1st.txt ${ImmPhen_DIR}/Data/Genotypes/checkQC/LabExMI_imputation_816x5699237_noMiss_noHet_noRel_noOutl_HWE_MAF_QCd ${ImmPhen_DataDir}/phenotypes_transformed.txt ${ImmPhen_DataDir}/covariates_selected.txt ${ImmPhen_DIR}/GWAS/regression ${PR_DIR}

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 10. Perform unconditional and conditional analyses on secondary trait (gene expression traits)
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Get list of genes per locus and create reference LD files
#SLURM
sbatch ${WORK_DIR}/Scripts/genes_and_LD.sh ${PR_DIR}/list.txt ${WORK_DIR} ${PR_DIR} ${WORK_DIR}/Scripts 
cat ${PR_DIR}/*/2nd/Pheno* > ${PR_DIR}/allgenes.txt
sort ${PR_DIR}/allgenes.txt | uniq > ${PR_DIR}/allgenes.txt_2
mv ${PR_DIR}/allgenes.txt_2 ${PR_DIR}/allgenes.txt

# Get 1MB window around genes and genes per cell type
Rscript ${WORK_DIR}/Scripts/get_gene_info.R ${PR_DIR}/allgenes.txt ${PR_DIR}/allgenes_info.txt
Rscript ${WORK_DIR}/Scripts/genes_per_celltype.R ${PR_DIR}/allgenes_info.txt ${WORK_DIR}/BP/Data/Pheno/phenotypes_transformed.txt ${PR_DIR}/allgenes_info_Tcells.txt
Rscript ${WORK_DIR}/Scripts/genes_per_celltype.R ${PR_DIR}/allgenes_info.txt ${WORK_DIR}/BP_Mono/Data/Pheno/phenotypes_transformed.txt ${PR_DIR}/allgenes_info_Monocytes.txt
Rscript ${WORK_DIR}/Scripts/genes_per_celltype.R ${PR_DIR}/allgenes_info.txt ${WORK_DIR}/BP_Neutro/Data/Pheno/phenotypes_transformed.txt ${PR_DIR}/allgenes_info_Neutrophils.txt

# Run conditional cis eQTL analyses on all genes
mkdir ${PR_DIR}/Tcells_assoc
mkdir ${PR_DIR}/Monocytes_assoc
mkdir ${PR_DIR}/Neutrophils_assoc
#!!! The sbatch --array information has to be updated in get_independent_2nd.sh for each cell types (wc -l ${PR_DIR}/allgenes_info_*.txt)
#SLURM
sbatch ${WORK_DIR}/Scripts/get_independent_2nd.sh ${PR_DIR}/allgenes_info_Tcells.txt ${WORK_DIR}/BP/Data/Geno/checkQC/all_noMiss_noHet_noRel_noOutl ${WORK_DIR}/BP/Data/Pheno/phenotypes_transformed.txt ${WORK_DIR}/BP/Data/Pheno/covariates_selected.txt ${WORK_DIR} ${PR_DIR}/Tcells_assoc 
#SLURM
sbatch ${WORK_DIR}/Scripts/get_independent_2nd.sh ${PR_DIR}/allgenes_info_Monocytes.txt ${WORK_DIR}/BP/Data/Geno/checkQC/all_noMiss_noHet_noRel_noOutl ${WORK_DIR}/BP_Mono/Data/Pheno/phenotypes_transformed.txt ${WORK_DIR}/BP_Mono/Data/Pheno/covariates_selected.txt ${WORK_DIR} ${PR_DIR}/Monocytes_assoc
#SLURM
sbatch ${WORK_DIR}/Scripts/get_independent_2nd.sh ${PR_DIR}/allgenes_info_Neutrophils.txt ${WORK_DIR}/BP/Data/Geno/checkQC/all_noMiss_noHet_noRel_noOutl ${WORK_DIR}/BP_Neutro/Data/Pheno/phenotypes_transformed.txt ${WORK_DIR}/BP_Neutro/Data/Pheno/covariates_selected.txt ${WORK_DIR} ${PR_DIR}/Neutrophils_assoc 

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 11. JLIM analysis
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Select all primary trait association statistics with a p < 1x10-5
bash ${WORK_DIR}/Scripts/select_assoc_1st.sh ${PR_DIR} 

# Select all secondary trait association statistics for each locus with p < 10-3
bash ${WORK_DIR}/Scripts/select_assoc_2nd.sh ${PR_DIR}

# Create association stat files and permutation files for secondary traits (gene expression traits) 
Rscript ${WORK_DIR}/Scripts/list_foruniqueperms.R ${PR_DIR}/list_1stand2ndtraits.txt ${PR_DIR}/list_1stand2ndtraits.txt_unique2nd
#SLURM
sbatch ${WORK_DIR}/Scripts/permutation2nd.sh ${PR_DIR} ${WORK_DIR}/BP/Data/Geno/checkQC/all_noMiss_noHet_noRel_noOutl ${WORK_DIR} 

# Run JLIM
#SLURM
sbatch ${WORK_DIR}/Scripts/Jlim.sh ${PR_DIR} ${WORK_DIR}/BP/Data/Geno/checkQC/all_noMiss_noHet_noRel_noOutl ${WORK_DIR} 
bash ${WORK_DIR}/Scripts/clean_JLIM.sh ${PR_DIR} ${WORK_DIR}   

# Create file missingJLIM.txt from list_1stand2ndtrait with lines that did not work and repeat those
# This is just done to be able to look in detail into the JLIM runs that "failed" - the ones for which no jlim.out.tsv file is created. They are rerun individually, to be able to look at the log files to find out, why the JLIM runs failed. Possible reasons for failure are 1) the minimal p-value for the secondary trait in the window JLIM selects automatically is too high (only SNPs shared by both data sets are included here) and 2) there are too few shared SNPs to run JLIM on this window. 
bash ${WORK_DIR}/Scripts/add_conditioninfo.sh ${PR_DIR} ${WORK_DIR}
Rscript ${WORK_DIR}/Scripts/get_missing_Jlim.R ${PR_DIR}/JLIMresults.txt ${PR_DIR}/list_1stand2ndtraits.txt ${PR_DIR}/missingJLIM.txt
#SLURM
sbatch ${WORK_DIR}/Scripts/Jlim_rep.sh ${PR_DIR} ${WORK_DIR}/BP/Data/Geno/checkQC/all_noMiss_noHet_noRel_noOutl ${WORK_DIR} 

# Remove failed JLIM run result files
bash ${WORK_DIR}/Scripts/clean_JLIM.sh ${PR_DIR} ${WORK_DIR}

# Combine JLIM results and add conditioning information
bash ${WORK_DIR}/Scripts/add_conditioninfo.sh ${PR_DIR} ${WORK_DIR}

# Info about which ones worked and failed (and why) can be found in ${PR_DIR}/info_JLIMS.txt. To create this the Log-files are evaluated. 
grep -l "is high" ${WORK_DIR}/Logs/PRS_runJLIM_21293627_* > ${PR_DIR}/high.txt
grep -l "few" ${WORK_DIR}/Logs/PRS_runJLIM_21293627_* > ${PR_DIR}/toofew.txt
Rscript ${WORK_DIR}/Scripts/info_missingJLIM.R ${PR_DIR}/list_1stand2ndtraits.txt ${PR_DIR}/JLIMresults.txt ${PR_DIR}/high.txt ${PR_DIR}/toofew.txt ${PR_DIR}/missingJLIM.txt 21293627 ${PR_DIR}/info_JLIMs.txt

# Select best JLIM results per ImmPhen/Gene combination
Rscript ${WORK_DIR}/Scripts/sort_JLIMresults.R ${PR_DIR}/JLIMresults.txt ${PR_DIR}/JLIMresults_selected.txt ${PR_DIR}/JLIMresults_selected_tested.txt

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 12. PRS analysis
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Prepare genetic data sets (get overlapping SNPs, remove ambiguous SNPs, flip SNPs when necessary)
bash ${WORK_DIR}/Scripts/prepare_data.sh ${WORK_DIR} ${PR_DIR}/PRS ${WORK_DIR}/BP/Data/Geno/checkQC/all_noMiss_noHet_noRel_noOutl ${ImmPhen_DIR}/Data/Genotypes/checkQC/LabExMI_imputation_816x5699237_noMiss_noHet_noRel_noOutl_HWE_MAF_QCd

# For each immune phenotype/gene expression trait combination select lead SNP to condition on
bash ${WORK_DIR}/Scripts/selectlead.sh ${PR_DIR}
awk '{ print $1"\t"$2"\t"$9 }' ${PR_DIR}/JLIMresults_selected_leads.txt >  ${PR_DIR}/JLIMresults_selected_leads_unique.txt
sort ${PR_DIR}/JLIMresults_selected_leads_unique.txt | uniq > ${PR_DIR}/JLIMresults_selected_leads_unique.txt_2
mv ${PR_DIR}/JLIMresults_selected_leads_unique.txt_2 ${PR_DIR}/JLIMresults_selected_leads_unique.txt

# Run genome-wide eQTL analysis  on genes conditioning on lead SNP
#SLURM
sbatch ${WORK_DIR}/Scripts/run_PRS_assoc_cond.sh ${WORK_DIR} ${PR_DIR} 

# Run unconditional genome-wide eQTL analysis on all genes
awk '{ print $1"\t"$2 }' ${PR_DIR}/JLIMresults_selected_leads.txt > ${PR_DIR}/JLIMresults_selected_genes_unique.txt
sort ${PR_DIR}/JLIMresults_selected_genes_unique.txt | uniq > ${PR_DIR}/JLIMresults_selected_genes_unique.txt_2
mv ${PR_DIR}/JLIMresults_selected_genes_unique.txt_2 ${PR_DIR}/JLIMresults_selected_genes_unique.txt
#SLURM
sbatch ${WORK_DIR}/Scripts/run_PRS_assoc.sh ${WORK_DIR} ${PR_DIR} 
 
# Calculate R2 for the immune phenotypes (using PRS calculated while conditioning on lead)
#SLURM
sbatch ${WORK_DIR}/Scripts/run_PRS_2pheno_condJLIM.sh ${PR_DIR}/JLIMresults_selected_complete.txt ${WORK_DIR} ${PR_DIR}

# Calculate R2 for the immune second phenotypes (using unconditional PRS)
#SLURM
sbatch ${WORK_DIR}/Scripts/run_PRS_2pheno_all.sh ${PR_DIR}/JLIMresults_selected_complete.txt

# Combine PRS results
#SLURM
sbatch ${WORK_DIR}/Scripts/combine_PRS.sh ${PR_DIR}/PRS condJLIM 
#SLURM
sbatch ${WORK_DIR}/Scripts/combine_PRS.sh ${PR_DIR}/MIP_BP/PRS all

# Print out results, preliminary plots 
for pthres in "5e-08" "1e-07" "5e-07" "1e-06" "5e-06" "1e-05" "5e-05" "0.0001" "0.001" "0.01"
do
Rscript ${WORK_DIR}/Scripts/plot_R2_PRS.R ${PR_DIR}/PRS/PRSresults_condJLIM_${pthres}.txt ${PR_DIR}/JLIMresults_selected_complete.txt condJLIM 0.00001 none ${PR_DIR}/PRS/plots_PRS_${pthres}_condJLIM.png ${PR_DIR}/PRS/PRS_res_${pthres}_condJLIM.txt
done
for pthres in "5e-08" "1e-07" "5e-07" "1e-06" "5e-06" "1e-05" "5e-05" "0.0001" "0.001" "0.01"
do
Rscript ${WORK_DIR}/Scripts/plot_R2_PRS.R ${PR_DIR}/PRS/PRSresults_all_${pthres}.txt ${PR_DIR}/JLIMresults_selected_complete.txt all 0.00001 none ${PR_DIR}/PRS/plots_PRS_${pthres}_all.png ${PR_DIR}/PRS/PRS_res_${pthres}_all.txt
done

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 13. TWMR analysis
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Get 1Mbp window for each gene
for DS in BP BP_Mono BP_Neutro
do
Rscript ${WORK_DIR}/Scripts/TWMR/define_windows.R ${WORK_DIR}/${DS}/Data/Pheno/PhenoList_all.txt ${WORK_DIR}/${DS}/Data/Pheno/Gtex/all.expression.bed ${WORK_DIR}/${DS}/Data/Pheno/PhenoList_all_windows.txt
done

# Copy genotype files that have been matched (shared SNPs between both data sets, flipped)
cp ${PR_DIR}/PRS/gen1.bim ${WORK_DIR}/BP/Data/Geno/checkQC/all_noMiss_noHet_noRel_noOutl_TWMR.bim
cp ${PR_DIR}/PRS/gen1.bed ${WORK_DIR}/BP/Data/Geno/checkQC/all_noMiss_noHet_noRel_noOutl_TWMR.bed
cp ${PR_DIR}/PRS/gen1.fam ${WORK_DIR}/BP/Data/Geno/checkQC/all_noMiss_noHet_noRel_noOutl_TWMR.fam 
cp ${PR_DIR}/PRS/gen2.bim ${WORK_DIR}/MIP/Data/Genotypes/checkQC/LabExMI_imputation_816x5699237_noMiss_noHet_noRel_noOutl_HWE_MAF_QCd_TWMR.bim
cp ${PR_DIR}/PRS/gen2.bed ${WORK_DIR}/MIP/Data/Genotypes/checkQC/LabExMI_imputation_816x5699237_noMiss_noHet_noRel_noOutl_HWE_MAF_QCd_TWMR.bed
cp ${PR_DIR}/PRS/gen2.fam ${WORK_DIR}/MIP/Data/Genotypes/checkQC/LabExMI_imputation_816x5699237_noMiss_noHet_noRel_noOutl_HWE_MAF_QCd_TWMR.fam

# Get SNPs in ImmPhen data set 
genfileImmPhen="${WORK_DIR}/MIP/Data/Genotypes/checkQC/LabExMI_imputation_816x5699237_noMiss_noHet_noRel_noOutl_HWE_MAF_QCd_TWMR"
awk '{print $2}' ${genfileImmPhen}.bim > ${genfileImmPhen}.snps

# Get SNPs with mismatch regarding A1 and A1
Rscript ${WORK_DIR}/Scripts/listSwitchedSNPs.R ${WORK_DIR}/MIP/Data/Genotypes/checkQC/LabExMI_imputation_816x5699237_noMiss_noHet_noRel_noOutl_HWE_MAF_QCd_TWMR.bim ${WORK_DIR}/BP/Data/Geno/checkQC/all_noMiss_noHet_noRel_noOutl_TWMR.bim ${PR_DIR}/MR/mismatch_snps.txt

# Run TWMR
pthres="1e-3"
#SLURM
sbatch ${WORK_DIR}/Scripts/TWMR/Ablauf.sh ${PR_DIR}/JLIMresults_selected.txt $pthres ${WORK_DIR} ${PR_DIR}/MR/TWMR ${genfileImmPhen}

# Combine results
bash ${WORK_DIR}/Scripts/TWMR/combine_results.sh ${PR_DIR}/MR/TWMR

# Plot results
Rscript ${WORK_DIR}/Scripts/TWMR/postTWMR.R ${PR_DIR}/MR/TWMR/all.alpha ${PR_DIR}/JLIMresults_selected_complete.txt ${PR_DIR}/MR/TWMR/TWMR_${pthres}.png


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
# 14. Two Sample MR 
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 

# Prepare ImmPhen assoc files 
awk ' (NR > 1) { if ($1 ~ /binary/) print $1; else print $1"_norm";}' ${PR_DIR}/JLIMresults_selected_complete.txt > ${PR_DIR}/MR/TwoSampleMR/ImmPhens.txt 
sort ${PR_DIR}/MR/TwoSampleMR/ImmPhens.txt | uniq > ${PR_DIR}/MR/TwoSampleMR/ImmPhens.txt_2 
mv ${PR_DIR}/MR/TwoSampleMR/ImmPhens.txt_2 ${PR_DIR}/MR/TwoSampleMR/ImmPhens.txt 
#SLURM
sbatch ${WORK_DIR}/Scripts/run_prepareMR_ImmPhens.sh ${PR_DIR}/MR/TwoSampleMR/ImmPhens.txt ${PR_DIR} ${WORK_DIR}

# Prepare Gene assoc files 
plink --bfile ${WORK_DIR}/MIP/Data/Genotypes/checkQC/LabExMI_imputation_816x5699237_noMiss_noHet_noRel_noOutl_HWE_MAF_QCd --freq --out ${WORK_DIR}/MIP/Data/Genotypes/checkQC/LabExMI_imputation_816x5699237_noMiss_noHet_noRel_noOutl_HWE_MAF_QCd
awk ' (NR > 1) {print $2"_"$3}' ${PR_DIR}/JLIMresults_selected_complete.txt > ${PR_DIR}/MR/TwoSampleMR/Genes.txt 
sort ${PR_DIR}/MR/TwoSampleMR/Genes.txt | uniq > ${PR_DIR}/MR/TwoSampleMR/Genes.txt_2 
mv ${PR_DIR}/MR/TwoSampleMR/Genes.txt_2 ${PR_DIR}/MR/TwoSampleMR/Genes.txt 
#SLURM
sbatch ${WORK_DIR}/Scripts/run_prepareMR_Genes.sh ${PR_DIR}/MR/TwoSampleMR/Genes.txt ${PR_DIR} ${WORK_DIR}

# Run MR
#SLURM
sbatch ${WORK_DIR}/Scripts/run_MR.sh ${PR_DIR}/JLIMresults_selected.txt ${PR_DIR}/PRS 

# Combine results
cat ${PR_DIR}/MR/TwoSampleMR/*results.txt > ${PR_DIR}/MR/TwoSampleMR/all.MRresults.txt
awk '{ if (NR == 1 || $1 != "id.exposure") print $0 }' ${PR_DIR}/MR/TwoSampleMR/all.MRresults.txt > ${PR_DIR}/MR/TwoSampleMR/all.MRresults.txt_2
mv ${PR_DIR}/MR/TwoSampleMR/all.MRresults.txt_2 ${PR_DIR}/MR/TwoSampleMR/all.MRresults.txt

# Plot Results
Rscript ${WORK_DIR}/Scripts/postTwoSampleMR.R ${PR_DIR}/MR/TwoSampleMR/all.MRresults.txt ${PR_DIR}/JLIMresults_selected_complete.txt "Inverse variance weighted" 0.00001 ${PR_DIR}/MR/TwoSampleMR/plots_IVW_0.00001.png 


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Some example Plots 
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# Example Plot for JLIM results

ImmPhen="N_CCR6pos_CD8pos.panel9"
Gene="ENSG00000012779.10"
DS="Neutrophils"
cond1st="uncond"
cond2nd="uncond"
chr=10
startbp=45778382
endbp=45978382
firsttraitp="Number of CCR6+CD8+ T cells"
firsttraitn="N_CCR6pos_CD8pos.panel9_norm"
phenoG="${WORK_DIR}/BP_Neutro/Data/Pheno/phenotypes_transformed.txt"

bash ${WORK_DIR}/Scripts/Submission/Figure_JLIM.sh $ImmPhen $Gene $DS $cond1st $cond2nd $chr $startbp $endbp "${firsttraitp}" $firsttraitn $phenoG ${PR_DIR} ${WORK_DIR} ${WORK_DIR}/MIP/GWAS/regression ${WORK_DIR}/MIP/Data/Genotypes/checkQC/LabExMI_imputation_816x5699237_noMiss_noHet_noRel_noOutl_HWE_MAF_QCd

# Example Plot for PRS results
bash ${WORK_DIR}/Scripts/Submission/Figure_PRS.sh ${WORK_DIR} ${PR_DIR}

# Example Plot for MR results
bash ${WORK_DIR}/Scripts/Submission/Figure_MR.sh ${WORK_DIR} ${PR_DIR}






