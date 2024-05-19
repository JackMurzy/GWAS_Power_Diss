# GWAS_Power_Diss
Code for my dissertation on quantifying the power of a GWAS.

## General Framework for Generating Genomic Data at Scale

The `generate_genomic_data` R function outputs genomic data for `i` individuals across `k` chromosomes with each containing `j` SNPs. From this data, `n` causal SNPs are assigned to generate `m` case phenotypes (the rest being controls). There are three key parts to the function:

1. Generating SNP genotype data
2. Assigning causal SNPs
3. Assigning case phenotypes







