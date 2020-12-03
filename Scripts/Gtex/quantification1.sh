#! /bin/bash


#SBATCH --job-name=run_pipeline
#SBATCH --output=/home/cg859/scratch60/Logs/STAR_run_%A_%a_log.txt
#SBATCH --partition=scavenge
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=10000
#SBATCH --time=10:00:00
#SBATCH --array=15-16

module load SAMtools

idsFile=$1
dir=$2
stardir=$3
wd=$4
filesFile=$5

#mkdir ${dir}/remove

readarray ids < $idsFile
readarrays files < $filesFile

index=$((SLURM_ARRAY_TASK_ID-1))

sample_id="${ids[$index]}"
sample_id="$(echo "$sample_id"|tr -d '\n')"
file_id="${files[$index]}"
file_id="$(echo "$file_id"|tr -d '\n')"

##ROSMAP:
#mv ${dir}/${sample_id}.r1.fastq.gz ${dir}/${sample_id}_1.fastq.gz
#mv ${dir}/${sample_id}.r2.fastq.gz ${dir}/${sample_id}_2.fastq.gz

##BP Tcells: 
#cp /gpfs/ysm/scratch60/kh593/blueprint/downloads/*_${sample_id}_* ${dir}/
#mv ${dir}/*${sample_id}*1.fastq.gz ${dir}/${sample_id}_1.fastq.gz
#mv ${dir}/*${sample_id}*2.fastq.gz ${dir}/${sample_id}_2.fastq.gz
#mkdir ${dir}/star_out

##BP Mono:
cp /gpfs/ysm/scratch60/kh593/blueprint/monocytes/completed/${file_id} 

#rm -f ${dir}/all.metrics.tsv > /dev/null


singularity exec ${wd}/Scripts/Gtex/gtex_pipeline.simg python3 /src/run_STAR.py ${dir}/${stardir} ${dir}/${sample_id}_1.fastq.gz ${dir}/${sample_id}_2.fastq.gz ${sample_id} --threads 10 --output_dir ${dir}/star_out

singularity exec ${wd}/Scripts/Gtex/gtex_pipeline.simg python3 /src/run_MarkDuplicates.py ${dir}/star_out/${sample_id}.Aligned.sortedByCoord.out.bam ${sample_id}.Aligned.sortedByCoord.out.md --output_dir ${dir} > ${dir}/${sample_id}.dups.log


if grep -q "Finished MarkDuplicates" ${dir}/${sample_id}.dups.log
then
	#rm ${dir}/${sample_id}_1.fastq.gz
	#rm ${dir}/${sample_id}_2.fastq.gz
	mv ${dir}/${sample_id}.dups.log ${dir}/Logs/
	rm -r ${dir}/star_out/${sample_id}*
	samtools index ${dir}/${sample_id}.Aligned.sortedByCoord.out.md.bam
	echo "${sample_id}" >> ${dir}/quant0done.txt
fi
















