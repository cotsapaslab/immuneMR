# immuneMR

## Contact
4 December 2020 <br>
[Christiane Gasperi](mailto:c.gasperi@tum.de) <br>
[Chris Cotsapas](mailto:cotsapas@broadinstitute.org) <br>

## Front matter
This README file documents the results of our study `Shared associations identify causal relationships between gene expression and immune cell phenotypes`(submitted). Here we provide an overview of the code provided in this distribution, the software dependencies, and the sources of genotype, gene expression and immune cell phenotype data required.   

## Disclaimers
* Throughout, we use the SLURM job scheduling system to distribute some of the more compute-intensive analyses on a cluster. We have tagged all instances with `#SLURM` in `MasterScript.sh` for ease of reference.
* The code here is provided as-is. We have validated that this runs without errors on our own systems, but do not guarantee that this will be the case in different computing environments.

## Data sets
We used gene expression data and genotype data from the BLUEPRINT consortium and immune cell phenotype and genotype data from the Milieu Intérieur project. 

### BLUEPRINT data set
+ Gene expression data: We downloaded FASTQ files for the RNAseq data for naive CD4+ T cells (169 individuals), CD14+ monocytes (193 individuals) and CD16+ neutrophils (196 individuals) from the BLUEPRINT consortium [Chen et al., Cell, 2016] (https://pubmed.ncbi.nlm.nih.gov/27863251). We used the GTEx pipeline for RNA-seq alignment, quantification and quality control (https://www.gtexportal.org/, Analysis Methods for V8). `Scripts/GenExp_phenodata.sh` gives an overview on how we applied the GTEx pipeline and can be used for orientation, but is not directly applicable to other data sets or computing environment needs user input and adaptations at various stages.
+ Genotype data: We obtained genotype data for all individuals for 7,008,524 variants acquired by whole genome sequencing (file name: EGAZ00001235598_release_vcf_06092016_All_chr.BPWP10_13_12_15.vcf). 

### Milieu Intérieur project data set
+ Phenotype data: We obtained flow cytometry measurements for 166 immune phenotypes of 816 individuals (file name: LabExMI_rawfacs.txt) and additional information including age, smoking status and CMV infection (file name: LabExMI_covariates.txt).
+ Genotype data: We obtained quality controlled and imputed genotype data from 816 individuals (file name: LabExMI_imputation_816x5699237.bim/.bed/.fam).

## Running our analysis pipeline
Here, we provide the code required to re-analyze the data and replicate our published results. For the seq alignment, quantification and quality control of the BLUEPRINT gene expression data please refer to the GTex pipeline (https://www.gtexportal.org/, Analysis Methods for V8).

### Pipeline overview
Our pipeline is split into fourteen distinct steps. These are:

+ Preparation of the Milieu Intérieur project phenotype data (quality control (QC) and normalization). 
+ QC of the Milieu Intérieur project genotype data
+ Preparation of the covariates for the Milieu Intérieur project data set
+ Alignment, quantification and quality control of the BLUEPRINT gene expression data
+ QC of the BLUEPRINT genotype data
+ Preparation of the covariates for the BLUEPRIN data set
+ GWAS analyses on the Milieu Intérieur project immune cell phenotypes 
+ Selection of genetic loci association with immune cell phenotypes
+ Perform conditional association analyses on the Milieu Intérieur project immune cell phenotypes 
+ Perform unconditional and conditional association analyses on the BLUEPRINT gene expression data (eQTL analysis)
+ JLIM analysis to identify shared association between gene expression and immune cell phenotypes 
+ Polygenic score (PGS) analysis to identify overall genetic correlation between gene expression and immune cell phenotypes 
+ Transcriptome-wide summary statistics-based Mendelian Randomization (TWMR) analysis (see Porcu et al., Nature Communications, 2019, https://www.nature.com/articles/s41467-019-10936-0)
+ Two sample MR (TSMR) analysis

### Software dependencies
+ PLINK v1.9
+ Pandoc 2.0.4
+ tabix 0.2.6
+ VCFtools 0.1.14
+ PRSice-2 (PRSice.R and PRSice_linux need to be put into Scripts/PRS) 
+ JLIM 2.0 (https://github.com/cotsapaslab/jlim, the jlim-master directory needs to be put into Scripts/JLIM)
+ R v3.4.1 or later (including packages ggplot2 v3.3.3, plyr v1.8.6, tidyr v1.0.2, gridExtra v.3, knitr v1.28, ggsci v2.9, gCMAP v1.22.1, readxl v1.3.1, png v0.1-7, biomaRt v2.34.2)
+ R v3.5.0 (used with the package TwosampleMR v0.4.25)
+ For software dependencies needed for the Gtex pipeline please see https://www.gtexportal.org/, Analysis Methods for V8.

### Running the pipeline
Except for the alignment, quantification and quality control of the BLUEPRINT gene expression data (for this see `Scripts/GenExp_phenodata.sh` and the Gtex pipeline) the entire pipeline can be run from `MasterScript.sh`, which calls the scripts included in this distribution to execute all fourteen steps of our pipeline. The script is not intended as a stand-alone, unsupervised piece of code - rather, it steps through the various procedures. User input is required at various stages.

### Description of files and directories in this distribution

| File/directory | Description |
| ------------- | ------------- |
| Scripts/  | Contains all scripts required in by the pipeline.  |
| MasterScript.sh  | Main script with fourteen QC and analysis steps of our pipeline.  |
| Scripts/GenExp_phenodata.sh  | Overview on how we applied the Gtex pipeline on the BLUEPRINT gene expression data.  |
| README.md  | This file.  |




