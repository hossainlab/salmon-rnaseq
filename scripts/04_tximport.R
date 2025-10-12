# RNA-seq Analysis in R: tximport and Gene-level Summarization
# Author: Md. Jubayer Hossain
# Affiliation: DeepBio Limited | CHIRAL Bangladesh
# Date: October 2025
# Description:
#   Imports transcript-level quantifications from Salmon
#   and summarizes to gene-level counts for DESeq2. 

# Install Bioconductor Packages 
BiocManager::install("tximport")
BiocManager::install("DESeq2")
BiocManager::install("EnsDb.Hsapiens.v86")

# Load libraries
library(tidyverse)
library(tximport)
library(DESeq2)
library(EnsDb.Hsapiens.v86)


# Get the mapping from transcript IDs to gene symbols 
# What are the columns in the database?
columns(EnsDb.Hsapiens.v86)

# Get the TXID and SYMBOL columns for all entries in database
tx2gene <- AnnotationDbi::select(EnsDb.Hsapiens.v86, 
                                   keys = keys(EnsDb.Hsapiens.v86),
                                   columns = c('TXID', 'SYMBOL'))

# Remove the gene ID column
tx2gene <- dplyr::select(tx2gene, -GENEID)

# Get the quant files and metadata
# Collect the sample quant files
samples <- list.dirs('outputs/Salmon_out', recursive = FALSE, full.names = FALSE)
samples

# check quant files 
quant_files <- file.path('outputs/Salmon_out', samples, 'quant.sf')
quant_files

# sample names 
names(quant_files) <- samples
print(quant_files)

# Ensure each file actually exists
# all should be TRUE
file.exists(quant_files)  

# Set up metadata frame
# Metadata for DESeq2: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778
col_data <- data.frame(
  row.names = samples,
  cell_line = rep(c("N61311","N052611","N080611","N061011"), each = 4),
  condition = rep(c("untreated","dexamethasone","albuterol","albuterol_dexamethasone"), times = 4)
)

# condition as factor 
col_data$condition <- factor(col_data$condition)


# Compile the tximport counts object and make DESeq dataset
# Get tximport counts object
txi <- tximport(files = quant_files, 
                type = 'salmon',
                tx2gene = tx2gene,
                ignoreTxVersion = TRUE)

# Make DESeq dataset
dds <- DESeqDataSetFromTximport(txi = txi,
                                colData = col_data,
                                design = ~condition)

# Do DESeq analysis
# PCA
vsd <- vst(dds)
plotPCA(vsd)

# DEG analysis
dds <- DESeq(dds)

# Get the results
resdf <- results(dds)

# MA plot 
plotMA(resdf)

# convert as data frame 
resdf <- as.data.frame(resdf)
resdf$gene <- rownames(resdf)
rownames(resdf) <- NULL

# Write to CSV file
write.csv(resdf, file = "outputs/DESeq2_results.csv", row.names = FALSE)

# Save as RDS (for reloading in R later)
saveRDS(resdf, file = "outputs/DESeq2_results.rds")
