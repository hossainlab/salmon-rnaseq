# load packages 
library(tidyverse)
library(RColorBrewer)
library(EnhancedVolcano)
library(pheatmap)
library(PoiClaClu)

# import DESeq2 results 
deseq2_results <- read.csv("outputs/DESeq2_results.csv")
names(deseq2_results)

# Volcano Plot
# Default cutoffs are log2FC > |2| and adjusted P-value < 0.05
EnhancedVolcano(
  deseq2_results, 
  lab = deseq2_results$gene, 
  x = "log2FoldChange", 
  y = "padj"
)

# Add custom log2FC and adjusted P-value cutoffs and size of points and labels 
EnhancedVolcano(
  deseq2_results, 
  lab = deseq2_results$gene, 
  x = "log2FoldChange", 
  y = "padj", 
  pCutoff = 0.001, 
  FCcutoff = 2, 
  pointSize = 1.5, 
  labSize = 3.0, 
  title = "Untreated vs. Treated"
)

# Adjust axis limits 
EnhancedVolcano(
  deseq2_results, 
  lab = deseq2_results$gene, 
  x = "log2FoldChange", 
  y = "padj", 
  pCutoff = 0.001, 
  FCcutoff = 2, 
  pointSize = 1.5, 
  labSize = 3.0, 
  xlim = c(-5, 5), 
  ylim = c(0, -log10(10e-10)), 
  title = "Untreated vs. Treated"
)

# Modify border and remove grid lines 
EnhancedVolcano(
  deseq2_results, 
  lab = deseq2_results$gene, 
  x = "log2FoldChange", 
  y = "padj", 
  pCutoff = 0.001, 
  FCcutoff = 2, 
  pointSize = 1.5, 
  labSize = 3.0, 
  xlim = c(-5, 5), 
  ylim = c(0, -log10(10e-10)), 
  border = "full", 
  borderWidth = 1.5, 
  borderColour = "black", 
  gridlines.major = FALSE, 
  title = "Untreated vs. Treated"
)