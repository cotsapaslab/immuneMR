#! /bin/bash


#SBATCH --job-name=Gtex_index
#SBATCH --output=/home/cg859/scratch60/Logs/STAR_index_%A_log.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --partition=scavenge
#SBATCH --mem-per-cpu=50000
#SBATCH --time=10:00:00

#module load RSEM

######Anpassen: 
# --sjbdOverhang = number of paired end basepair reads -1:
# für CMC: "One-hundred base pair paired end reads were obtained on a HiSeq 2500." --> num=100 (99??)
# für ROSMAP: "Sequencing was carried out using the Illumina HiSeq2000 with 101 bp paired end reads for a targeted coverage of 50M paired reads. --> num=100



dir=$1
num=$2
wd=$3

mkdir ${dir}/star_index_overhang${num}



singularity exec ${wd}/Scripts/Gtex/gtex_pipeline.simg STAR --runMode genomeGenerate --genomeDir ${dir}/star_index_overhang${num} --genomeFastaFiles ${dir}/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta --sjdbGTFfile ${dir}/gencode.v26.GRCh38.annotation.ERCC.gtf --sjdbOverhang ${num} --runThreadN 10

#singularity exec /home/cg859/scratch60/GtexPipeline/gtex_pipeline.simg rsem-prepare-reference --gtf /home/cg859/scratch60/GtexPipeline/references/gencode.v26.GRCh38.annotation.ERCC.gtf --num-threads 4 /home/cg859/scratch60/GtexPipeline/references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta /home/cg859/scratch60/GtexPipeline/rsem_reference/rsem_reference




