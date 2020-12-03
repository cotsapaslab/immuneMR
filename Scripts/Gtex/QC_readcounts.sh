#! /bin/bash

#SBATCH --job-name=QC_Gtex
#SBATCH --output=/home/cg859/scratch60/Logs/STAR_QC_%A_log.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --partition=scavenge
#SBATCH --mem-per-cpu=10000
#SBATCH --time=10:00:00


module load PLINK/1.90-beta4.6
module load SAMtools
module load tabix

dir=$1
wd=$2


# Normalization
singularity exec ${wd}/Scripts/Gtex/gtex_pipeline_2.simg python3 /src/eqtl_prepare_expression.py ${dir}/all.gene_tpm_QC1.gct.gz ${dir}/all.gene_reads_QC1.gct.gz ${dir}/gencode.v26.GRCh38.ERCC.genes.gtf ${dir}/sampleLookup.txt ${dir}/all_chroms.list ${dir}/all --tpm_threshold 0.1 --count_threshold 6 --sample_frac_threshold 0.2 --normalization_method tmm

# PEER calculation
cd ${dir}
singularity exec ${wd}/Scripts/Gtex/gtex_pipeline_2.simg Rscript /src/run_PEER.R ${dir}/all.expression.bed.gz all 60




