#!/bin/sh

gunzip *.gz
cat $(find . -mindepth 1 -maxdepth 1 -name \*.logistic | sort -V) >all.plink
grep -v SNP all.plink >all2.plink
fn=$(ls *.assoc.logistic | head -1)
head -1 ${fn} >header.txt
cat header.txt all2.plink >all3.plink
mv all3.plink ${1}.logistic
rm all2.plink
rm all.plink
mv ${1}.logistic ..
gzip ../${1}.logistic
fn=$(ls *.log | head -1)
cp ${fn} ../${1}.log
gzip *.logistic
rm header.txt
