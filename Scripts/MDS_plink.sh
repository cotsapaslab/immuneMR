#! /bin/bash


rem=$1
file=$2
prefix=$3

module load PLINK/1.90-beta4.6

plink --bfile ${file}/${prefix} --exclude range ${rem} --geno 0.02 --hwe 1e-3 --maf 0.05 --set-hh-missing --indep-pairwise 200 100 0.2 --out ${file}/${prefix}_prunelist
plink --bfile ${file}/${prefix} --extract ${file}/${prefix}_prunelist.prune.in --cluster --mds-plot 10 eigendecomp --out ${file}/${prefix}_cluster
plink --bfile ${file}/${prefix} --extract ${file}/${prefix}_prunelist.prune.in --cluster --pca 10 header --out ${file}/${prefix}_cluster

rm ${file}/${prefix}_prunelist*
rm ${file}/${prefix}_cluster.cluster*




