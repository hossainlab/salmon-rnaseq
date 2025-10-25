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


# Get the mapping from transcript IDs to gene symbols 
# What are the columns in the database?
columns(EnsDb.Hsapiens.v86)

# Get the TXID and SYMBOL columns for all entries in database
tx2gene <- AnnotationDbi::select(EnsDb.Hsapiens.v86, 
                                 keys = keys(EnsDb.Hsapiens.v86),
                                 columns = c('TXID', 'SYMBOL'))

# Remove the gene ID column
tx2gene <- dplyr::select(tx2gene, -GENEID)


# Compile the tximport counts object and make DESeq dataset
# Get tximport counts object
txi <- tximport(files = quant_files, 
                type = 'salmon',
                tx2gene = tx2gene,
                ignoreTxVersion = TRUE)

# class of txi 
class(txi)

# raw counts 
raw_counts <- txi$counts
write.csv(raw_counts, "outputs/tables/raw_counts.csv", row.names = FALSE)
write_rds(raw_counts, "outputs/tables/raw_counts.rds")


# TPM 
tpm_counts <- txi$abundance
write.csv(tpm_counts, "outputs/tables/tpm_counts.csv", row.names = FALSE)
write_rds(tpm_counts, "outputs/tables/tpm_counts.rds")


# Make DESeq dataset
dds <- DESeqDataSetFromTximport(txi = txi,
                                colData = col_data,
                                design = ~condition)

# Principal Component Analysis 
rlog_dds <- rlog(dds)

pca_data <- plotPCA(rlog_dds, intgroup = "condition", returnData = TRUE)
write_rds(pca_data, "outputs/tables/pca_data.rds")

# Differential Gene Expression Analysis 
dds <- DESeq(dds)

# Get the results
resdf <- results(dds)

write_rds(resdf, "outputs/tables/res_dds.rds")


# MA plot 
plotMA(resdf)

# convert as data frame 
resdf <- as.data.frame(resdf)

resdf$gene <- rownames(resdf)
rownames(resdf) <- NULL

# Save as RDS (for reloading in R later)
saveRDS(resdf, "outputs/tables/DESeq2_results.rds")