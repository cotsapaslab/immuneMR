# immuneMR

## Contact
4 Decmber 2020 <br>
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


