#! /bin/bash


#SBATCH --job-name=run_pipeline
#SBATCH --output=/home/cg859/scratch60/Logs/STAR_combine_%A_log.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --partition=scavenge
#SBATCH --mem-per-cpu=10000
#SBATCH --time=10:00:00


idsFile=$1
dir=$2
wd=$3
num=$4

readarray ids < $idsFile
numids=${#ids[@]}

listreads=""
listrpkm=""
listexons=""

for (( i=0; i<$numids; i++))
do

sample_id="${ids[$i]}"
sample_id="$(echo "$sample_id"|tr -d '\n')"

reads="${dir}/${sample_id}.gene_reads.gct.gz"
rpkm="${dir}/${sample_id}.gene_rpkm.gct.gz"
exons="${dir}/${sample_id}.exon_reads.gct.gz"

listreads="${listreads} ${reads}"
listrpkm="${listrpkm} ${rpkm}"
listexons="${listexons} ${exons}"
done


echo $listreads

singularity exec ${wd}/Scripts/Gtex/gtex_pipeline.simg python3 /src/combine_GCTs.py ${listreads} ${dir}/all${num}.gene_reads
singularity exec ${wd}/Scripts/Gtex/gtex_pipeline.simg python3 /src/combine_GCTs.py ${listrpkm} ${dir}/all${num}.gene_rpkm 
singularity exec ${wd}/Scripts/Gtex/gtex_pipeline.simg python3 /src/combine_GCTs.py ${listexons} ${dir}/all${num}.exon_rpkm 












