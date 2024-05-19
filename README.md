# GWAS_Power_Diss
This is the repository for the code behind my dissertation on Quantifying the Power of GWAS. It includes code to simulate genetic data in a fast and simple manner. As well as code for interpreting results and generating HTML reports.

## Files
- [genetic_sim.R](./genetic_sim.R): This file contains the main simulation logic. It is split into three sub-functions designed to simulate genomic data, assign causal SNPs and assign phenotypes
- [genetic_assoc_stats_function.R](./genetic_assoc_stats_function.R): This file contains a plethora of functions designed to generate association statistics for a PLINK --assoc GWAS output. These include power and inflation calculations, Bonferroni and FDR thresholds, F-values and ROC curves.
- [gwas_html_boiler.Rmd](./gwas_html_boiler.Rmd): This is the markdown template for the per-run HTML reports generated. This tracks stats such as run parameters, manhattan and qq plots as well as ROC curves and general association stats. 
- [summary_gwas_html_boiler.Rmd](./summary_gwas_html_boiler.Rmd): This is the markdown template for the per-test summary HTML reports. This plots power, F-scores, AUC and inflation statistics across ALL runs. 
- [ggqqman.R](./ggqqman.R): A custom function to generate qq and manhattan plots using ggplot2.
- [gen_ld_matrix.R](./gen_ld_matrix.R): A custom function to calculate the LD matrix between SNPs based off genotype data. This uses the statistical approach outlined by [T.J.Hui and A.Burt, 2020](https://bmcgenomdata.biomedcentral.com/articles/10.1186/s12863-020-0818-9#citeas). 


## General Framework for Generating Genomic Data at Scale

The `genetic_sim` R function outputs genomic data for `i` individuals across `k` chromosomes with each containing `j` SNPs. From this data, `n` causal SNPs are assigned to generate `m` case phenotypes (the rest being controls). There are three key parts to the function:

1. Generating SNP genotype data
2. Assigning causal SNPs
3. Assigning case phenotypes







