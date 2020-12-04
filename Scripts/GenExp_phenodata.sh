

############################################################
#We downloaded FASTQ files and used the GTEx pipeline for RNA-seq alignment, quantification and quality control (see and https://www.gtexportal.org/, Analysis Methods for V8 and https://github.com/broadinstitute/gtex-pipeline). 
#This script is intended to give an overview on how we applied the GTEx pipeline and can be used for orientation, but is not directly applicable to other data sets and needs user input at various stages. Here, the pipeline is used on the BluePrint Tcell data. 
############################################################
 
###########################################################
# Steps performed to set up the environment to run the Gtex pipeline: 
#1) Download gencode.v26.annotation.gtf from ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/
#2) python3 collapse_annotation.py gencode.v26.annotation.gtf gencode.v26.GRCh38.genes.gtf  (Script can be found at https://github.com/broadinstitute/gtex-pipeline/blob/master/gene_model/collapse_annotation.py)
#3) mv gencode.v26.annotation.gtf gencode.v26.GRCh38.annotation.gtf
#4) Download ERCC spike-in reference annotations from:  https://tools.thermofisher.com/content/sfs/manuals/ERCC92.zip
#5) sed 's/ERCC-/ERCC_/g' ERCC92.fa > ERCC92.patched.fai
#6) python3 changeERCC.py (created from https://github.com/broadinstitute/gtex-pipeline/blob/master/TOPMed_RNAseq_pipeline.md))
#7) cat gencode.v26.GRCh38.annotation.gtf ERCC92.genes.patched.gtf > gencode.v26.GRCh38.annotation.ERCC.gtf
#8) cat gencode.v26.GRCh38.genes.gtf ERCC92.genes.patched.gtf > gencode.v26.GRCh38.ERCC.genes.gtf
#9) Download Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.tar.gz from https://github.com/broadinstitute/gtex-pipeline/blob/master/TOPMed_RNAseq_pipeline.md
#10) mkdir ${WORK_DIR}/Gtexreferences
#11) mv gencode* Gtexreferences/ \ mv Homo_sapiens* Gtexreferences/ \ mv ERCC* Gtexreferences/
#12) https://github.com/francois-a/rnaseqc/releases/download/v1.1.9/RNA-SeQC_1.1.9.zip
#13) module load STAR; module load RSEM; module load SAMtools
#14) singularity build gtex_pipeline_2.simg docker://broadinstitute/gtex_eqtl:V8
#15) singularity build gtex_pipeline.simg docker://broadinstitute/gtex_rnaseq:V8
#16) sbatch Scripts/STARindex_RSEMref.sh #manually change parameters in script
#17) sbatch Scripts/run_pipeline.sh #manually change parameters in script
###########################################################


module load PLINK/1.90-beta4.6
module load R
module load Pandoc
module load VCFtools
module load tabix

WORK_DIR="~/home"


DS="BP"
DATA_DIR="${WORK_DIR}/${DS}/Data"

mkdir ${WORK_DIR}/${DS}/Data/Pheno/Gtex/
	 
# copy reference files for Gtex-Pipeline
cp ${WORK_DIR}/Gtexreferences/* ${WORK_DIR}/${DS}/Data/Pheno/Gtex/

# prepare STAR index
#SLURM
sbatch ${WORK_DIR}/Scripts/Gtex/STARindex_RSEMref.sh ${WORK_DIR}/${DS}/Data/Pheno/Gtex 100 ${WORK_DIR}

###### copy FASTQ-Files from BluePrint and create IDlist.txt (for BP_Mono and BP_Neutro use quant1_single.....sh and create two lists File_List_ID.txt and File_List_File.txt

# run GTex Pipeline part 1
#SLURM
sbatch ${WORK_DIR}/Scripts/Gtex/quantification1.sh ${WORK_DIR}/${DS}/Data/Pheno/Gtex/IDlist.txt ${WORK_DIR}/${DS}/Data/Pheno/Gtex star_index_overhang99 ${WORK_DIR}
#SLURM
sbatch ${WORK_DIR}/Scripts/Gtex/quantification2.sh ${WORK_DIR}/${DS}/Data/Pheno/Gtex/IDlist_part1_use.txt ${WORK_DIR}/${DS}/Data/Pheno/Gtex star_index_overhang99 ${WORK_DIR}
mv ${WORK_DIR}/Logs/STAR* ${WORK_DIR}/Logs/${DS}/


# Merge single sample output
#SLURM
sbatch ${WORK_DIR}/Scripts/Gtex/combine_output.sh ${WORK_DIR}/${DS}/Data/Pheno/Gtex/IDlist.txt ${WORK_DIR}/${DS}/Data/Pheno/Gtex/alldone ${WORK_DIR} ${nummm}

# Combine Metrics-files
bash ${WORK_DIR}/Scripts/Gtex/combine_metrics.sh ${DATA_DIR}/Pheno/Gtex/IDlist.txt ${DATA_DIR}/Pheno/Gtex

# Create TPM file
gunzip ${DATA_DIR}/Pheno/Gtex/all.gene_rpkm.gct.gz
Rscript ${WORK_DIR}/Scripts/Gtex/rpkm_to_tpm.R ${DATA_DIR}/Pheno/Gtex/all.gene_rpkm.gct ${DATA_DIR}/Pheno/Gtex/all.gene_tpm.gct
head -n 2 ${DATA_DIR}/Pheno/Gtex/all.gene_rpkm.gct > ${DATA_DIR}/Pheno/Gtex/header.gct
cat ${DATA_DIR}/Pheno/Gtex/header.gct ${DATA_DIR}/Pheno/Gtex/all.gene_tpm.gct > ${DATA_DIR}/Pheno/Gtex/all.gene_tpm.gct_2
mv ${DATA_DIR}/Pheno/Gtex/all.gene_tpm.gct_2 ${DATA_DIR}/Pheno/Gtex/all.gene_tpm.gct
rm ${DATA_DIR}/Pheno/Gtex/header.gct
	
#Sample QC
# Remove samples with less than 10000000 reads 
gunzip ${DATA_DIR}/Pheno/Gtex/all.gene_reads.gct.gz
head -n 2 ${DATA_DIR}/Pheno/Gtex/all.gene_rpkm.gct > ${DATA_DIR}/Pheno/Gtex/rpkm.header
head -n 2 ${DATA_DIR}/Pheno/Gtex/all.gene_reads.gct > ${DATA_DIR}/Pheno/Gtex/reads.header
head -n 2 ${DATA_DIR}/Pheno/Gtex/all.gene_tpm.gct > ${DATA_DIR}/Pheno/Gtex/tpm.header
Rscript ${WORK_DIR}/Scripts/Gtex/removeSamplesReads.R ${DATA_DIR}/Pheno/Gtex/all.metrics.tsv ${DATA_DIR}/Pheno/Gtex/all.gene_reads.gct ${DATA_DIR}/Pheno/Gtex/all.gene_rpkm.gct ${DATA_DIR}/Pheno/Gtex/all.gene_tpm.gct
cat ${DATA_DIR}/Pheno/Gtex/rpkm.header ${DATA_DIR}/Pheno/Gtex/all.gene_rpkm.gct_QC_tmp > ${DATA_DIR}/Pheno/Gtex/all.gene_rpkm_QC1.gct
cat ${DATA_DIR}/Pheno/Gtex/reads.header ${DATA_DIR}/Pheno/Gtex/all.gene_reads.gct_QC_tmp > ${DATA_DIR}/Pheno/Gtex/all.gene_reads_QC1.gct
cat ${DATA_DIR}/Pheno/Gtex/tpm.header ${DATA_DIR}/Pheno/Gtex/all.gene_tpm.gct_QC_tmp >  ${DATA_DIR}/Pheno/Gtex/all.gene_tpm_QC1.gct
gzip ${DATA_DIR}/Pheno/Gtex/all.*.gct
rm ${DATA_DIR}/Pheno/Gtex/*tmp
rm ${DATA_DIR}/Pheno/Gtex/*.header

#create chrom_list File with chr1\nchr2\nchr3...
paste -d"\t" ${DATA_DIR}/Pheno/Gtex/IDlist.txt ${DATA_DIR}/Pheno/Gtex/IDlist.txt > ${DATA_DIR}/Pheno/Gtex/sampleLookup.txt
printf "sample_id\tparticipant_id\n" > ${DATA_DIR}/Pheno/Gtex/sampleLookup.txt_header
cat ${DATA_DIR}/Pheno/Gtex/sampleLookup.txt_header ${DATA_DIR}/Pheno/Gtex/sampleLookup.txt > ${DATA_DIR}/Pheno/Gtex/sampleLookup.txt_2
mv ${DATA_DIR}/Pheno/Gtex/sampleLookup.txt_2 ${DATA_DIR}/Pheno/Gtex/sampleLookup.txt
rm ${DATA_DIR}/Pheno/Gtex/sampleLookup.txt_header
#SLURM
sbatch ${WORK_DIR}/Scripts/Gtex/QC_readcounts.sh ${DATA_DIR}/Pheno/Gtex ${WORK_DIR}

# Bed to Phenotype-Format
head -n 1 ${DATA_DIR}/Pheno/Gtex/all.expression.bed > ${DATA_DIR}/Pheno/Gtex/samples.txt
tr '\t' '\n' <${DATA_DIR}/Pheno/Gtex/samples.txt > ${DATA_DIR}/Pheno/Gtex/samples2.txt
mv ${DATA_DIR}/Pheno/Gtex/samples2.txt ${DATA_DIR}/Pheno/Gtex/samples.txt
Rscript ${WORK_DIR}/Scripts/Gtex/BedtoPheno.R ${DATA_DIR}/Pheno/Gtex/all.expression.bed ${DATA_DIR}/Pheno/Gtex/samples.txt ${DATA_DIR}/Geno/checkQC/all_noMiss_noHet_noRel_noOutl.fam ${DATA_DIR}/Pheno/phenotypes_transformed.txt CMC

# Select list of genes per focus
Rscript ${WORK_DIR}/Scripts/phenodata/genesperlocus.R ${MSGWAS_DIR}/MS-indexSNP.tsv ${DATA_DIR}/Pheno/Gtex/all.expression.bed ${DATA_DIR}/Pheno/

