#! /bin/bash


#SBATCH --job-name=run_pipeline
#SBATCH --output=/home/cg859/scratch60/Logs/STAR_run_%A_%a_log.txt
#SBATCH --cpus-per-task=5
#SBATCH --partition=general
#SBATCH --mem-per-cpu=20000
#SBATCH --time=12:00:00
#SBATCH --array=1-109


module load SAMtools

idsFile=$1
dir=$2
stardir=$3
wd=$4



readarray ids < $idsFile

index=$((SLURM_ARRAY_TASK_ID-1))

sample_id="${ids[$index]}"
sample_id="$(echo "$sample_id"|tr -d '\n')"

mkdir ${dir}/${sample_id}_run

cp ${dir}/gencode.v26.GRCh38.ERCC.genes.gtf ${dir}/${sample_id}_run/.
cp ${dir}/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta ${dir}/${sample_id}_run/.
cp ${dir}/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta.fai ${dir}/${sample_id}_run/.
mv ${dir}/${sample_id}.Aligned.sortedByCoord.out.md.bam ${dir}/${sample_id}_run/
mv ${dir}/${sample_id}.Aligned.sortedByCoord.out.md.bam.bai ${dir}/${sample_id}_run/

singularity exec ${wd}/Scripts/Gtex/gtex_pipeline.simg /src/run_rnaseqc.py ${sample_id}.Aligned.sortedByCoord.out.md.bam gencode.v26.GRCh38.ERCC.genes.gtf Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta ${sample_id} --output_dir ${dir}/${sample_id}_run -m 500 --java /usr/lib/jvm/java-1.7.0-openjdk-amd64/bin/java > ${dir}/${sample_id}.rnaseqc.log
rm ${dir}/${sample_id}_run/gencode.v26.GRCh38.ERCC.genes.gtf 
rm ${dir}/${sample_id}_run/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta
mv ${dir}/${sample_id}_run/* ${dir}/
rm -r ${dir}/${sample_id}_run

if grep -q "Finished Successfully" ${dir}/${sample_id}.rnaseqc.log
then
	if [ -f ${dir}/${sample_id}.exon_reads.gct.gz ]
	then
		if grep -q "Finished RNA-SeQC" ${dir}/${sample_id}.rnaseqc.log
		then
		rm ${dir}/${sample_id}.Alig*
		mv ${dir}/${sample_id}.* ${dir}/alldone/ 
		mv ${dir}/${sample_id}/ ${dir}/alldone/
		echo "${sample_id}" >> ${dir}/quant1done.txt
		fi
	fi
fi





















