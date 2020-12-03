#! /bin/bash


WORK_DIR=$1
PR_DIR=$2

module load R

Rscript ${WORK_DIR}/Scripts/Submission/plot_Figure_MR.R ${PR_DIR}/MR/TWMR/all.alpha ${PR_DIR}/MR/TwoSampleMR/all.MRresults.txt ${PR_DIR}/JLIMresults_selected_complete.txt "Inverse variance weighted" 0.00001 ${WORK_DIR}/Figure_MR.png

