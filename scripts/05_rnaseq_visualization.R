# load packages 
library(tidyverse)
library(tidyplots)
library(RColorBrewer)
library(EnhancedVolcano)
library(pheatmap)
library(PoiClaClu)

# import pca data 
pca_data <- readRDS("outputs/tables/pca_data.rds")
percentVar <- round(100 * attr(pca_data, "percentVar"), 1)

pca_data |>
  tidyplot(x = PC1, y = PC2, color = group) |>
  add_data_points(size = 1.8, white_border = TRUE) |>
  add_ellipse(level = 0.68, linewidth = 0.6) |>
  adjust_x_axis_title(paste0("PC1 (", percentVar[1], "%)")) |>
  adjust_y_axis_title(paste0("PC2 (", percentVar[2], "%)")) |>
  adjust_colors(colors_discrete_apple) |>
  adjust_legend_title("Group")
ggsave("outputs/figures/PCA_Plot.jpg", width = 7, height = 5.5, dpi = 300)
ggsave("outputs/figures/PCA_Plot.pdf")


# Volcano plot 
df <- 
  readRDS("outputs/tables/DESeq2_results.rds") |> 
  mutate(
    neg_log10_padj = -log10(padj),
    direction = if_else(log2FoldChange > 0, "up", "down", NA),
    candidate = abs(log2FoldChange) >= 1 & padj < 0.05
  )

df |> 
  tidyplot(x = log2FoldChange, y = neg_log10_padj) |> 
  add_data_points(data = filter_rows(!candidate),
                  color = "lightgrey", rasterize = TRUE) |> 
  add_data_points(data = filter_rows(candidate, direction == "up"),
                  color = "#FF7777", alpha = 0.5) |> 
  add_data_points(data = filter_rows(candidate, direction == "down"),
                  color = "#7DA8E6", alpha = 0.5) |> 
  add_reference_lines(x = c(-1, 1), y = -log10(0.05)) |> 
  add_data_labels_repel(data = min_rows(padj, 6, by = direction), label = gene,
                        color = "#000000", min.segment.length = 0, background = TRUE) |> 
  adjust_x_axis_title("$Log[2]~fold~change$") |> 
  adjust_y_axis_title("$-Log[10]~italic(P)~adjusted$")

ggsave("outputs/figures/Volcano_Plot.pdf")


# EnhancedVolcano
EnhancedVolcano(
  df, 
  lab = df$gene, 
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
ggsave("outputs/figures/EnhancedVolcano_Plot.pdf")




