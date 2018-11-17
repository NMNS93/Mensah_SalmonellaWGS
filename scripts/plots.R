# plots.R
# Generate figures for manuscript. 
#     Inputs: None
#     Outputs: 
#         plots/gene_freqs.png: A bar graph of gene frequencies

# Import functions and libraries
source('scripts/functions.R')

# Import plot data
raw_genotype_freqs = read.csv('results/tables/genotype_counts.csv',colClasses=c("NULL",NA,NA))