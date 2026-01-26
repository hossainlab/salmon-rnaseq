# Install package
pak::pkg_install(c(
  "EnhancedVolcano", 
  "pheatmap",
  "PoiClaClu",
  "tidyplots",
  "visdat", 
  "naniar"
))

# Load packages 
library(tidyverse)
library(tidyplots)
library(RColorBrewer)


# read pca data 
pca_data <- read_rds("outputs/tables/pca_data.rds")

# principle component analysis 
pca_data |> 
  tidyplot(x = PC1, y = PC2, color = group) |> 
  add_data_points(size = 1.3, white_border = TRUE) |> 
  add_ellipse() |> 
  adjust_x_axis_title(paste0("Component 1 (", round(pca_data$pc1_var*100, digits = 1), "%)")) |> 
  adjust_y_axis_title(paste0("Component 2 (", round(pca_data$pc2_var*100, digits = 1), "%)")) |> 
  adjust_colors(colors_discrete_apple) |> 
  adjust_legend_title("Group")
ggsave("outputs/figures/PCA_Plot.pdf")

# read deseq2 results 
deseq2_results <- read_rds("outputs/tables/DESeq2_results.rds")

# check missing value 
sum(is.na(deseq2_results))

# visualize missing values 
visdat::vis_dat(deseq2_results)
visdat::vis_miss(deseq2_results)

# which vars
naniar::miss_var_which(deseq2_results)
naniar::gg_miss_which(deseq2_results)

# remove missing values 
deseq2_results_clean <- deseq2_results |> 
  drop_na()

# export clean data 
write_rds(deseq2_results_clean, "outputs/tables/deseq2_results_clean.rds")

# read clean data 
deseq2_clean <- read_rds("outputs/tables/deseq2_results_clean.rds")

# filter significant genes 
sigs_genes <- deseq2_clean |> 
  dplyr::filter(padj < 0.05)

# filter up-regulated 
up_genes <- deseq2_clean |> 
  dplyr::filter(padj < 0.05 & log2FoldChange > 1)

# filter low expressed up-regulated genes 
up_genes |> 
  dplyr::filter(baseMean > 10)

# filter down-regulated 
down_genes <- deseq2_clean |> 
  dplyr::filter(padj < 0.05 & log2FoldChange < -1)

# filter low expressed down-regulated genes 
down_genes |> 
  dplyr::filter(baseMean > 10)

# combine into a summary 
summary_degs <- deseq2_clean |> 
  dplyr::mutate(
    regulation = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Up", 
      padj < 0.05 & log2FoldChange < -1 ~ "Down", 
      TRUE ~ "NS" # not significant 
    )
)



