# load packages 
library(tidyverse)
library(RColorBrewer)
library(EnhancedVolcano)
library(pheatmap)
library(PoiClaClu)

# import DESeq2 results 
deseq2_results <- read.csv("outputs/DESeq2_results.csv")
names(deseq2_results)

# Remove any genes that do not have any entrez identifiers
results_sig_entrez <- subset(deseq2_results, is.na(entrez) == FALSE)

# Create a matrix of gene log2 fold changes
gene_matrix <- results_sig_entrez$log2FoldChange

# Add the entrezID's as names for each logFC entry
names(gene_matrix) <- results_sig_entrez$entrez

# View the format of the gene matrix
##- Names = ENTREZ ID
##- Values = Log2 Fold changes
head(gene_matrix)