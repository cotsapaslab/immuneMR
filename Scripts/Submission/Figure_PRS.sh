#! /bin/bash


module load R

WORK_DIR=$1
PR_DIR=$2

pthres1="5e-08"
pthres2="5e-05"
Rscript ${WORK_DIR}/Scripts/Submission/plot_Figure_PRS.R ${PR_DIR}/PRS/PRSresults_condJLIM_${pthres1}.txt ${PR_DIR}/JLIMresults_selected_complete.txt condJLIM 0.00001 ${PR_DIR}/PRS/PRSresults_condJLIM_${pthres2}.txt ${PR_DIR}/JLIMresults_selected_complete.txt condJLIM 0.00001 ${WORK_DIR}/Figure_PRS.png




