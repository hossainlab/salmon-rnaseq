# Install clusterProfiler
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("enrichplot")
BiocManager::install("GOSemSim")


# load packages 
library(tidyverse)
library(tidyplots)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# Finding Pathways from Differential Expressed Genes
degs_df <- readRDS("outputs/tables/DESeq2_results.rds")

# Background universe = all tested genes (after any prefiltering you used before DESeq)
universe_symbols <- degs_df  |> 
  filter(!is.na(gene))  |> 
  pull(gene)  |> 
  unique()

# Significant DE genes (tune thresholds as needed)
sig_genes <- degs_df  |> 
  filter(!is.na(padj))  |> 
  filter(padj <= 0.05, abs(log2FoldChange) >= 0.58)  |>   # |LFC| >= 0.58 ≈ 1.5x
  pull(gene)  |> 
  unique()

# Map sig and universe to ENTREZ IDs
map_sig <- bitr(sig_genes, 
                fromType="SYMBOL", 
                toType="ENTREZID", 
                OrgDb=org.Hs.eg.db)


map_univ <- bitr(universe_symbols, 
                 fromType="SYMBOL", 
                 toType="ENTREZID", 
                 OrgDb=org.Hs.eg.db)

sig_entrez <- unique(map_sig$ENTREZID)
univ_entrez <- unique(map_univ$ENTREZID)

# quick sanity checks
length(sig_entrez)
length(univ_entrez)  


# KEGG Over-Representation Analysis (ORA)
ekegg <- enrichKEGG(
  gene          = sig_entrez,
  universe      = univ_entrez,
  organism      = "hsa",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.10,
  pAdjustMethod = "BH",
  minGSSize     = 10,
  maxGSSize     = 500
)

# View & save
head(as.data.frame(ekegg), 10)
write_csv(as.data.frame(ekegg), "outputs/tables/kegg_ORA_results.csv")

# Plot (publication-ready)
p_dot <- dotplot(ekegg, showCategory = 20) +
  ggtitle("KEGG pathway enrichment (ORA)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face="bold", hjust=0.5))
ggsave("outputs/figures/kegg_ORA_dotplot.png", p_dot, width=7, height=5, dpi=300)






