rm(list=ls())
setwd("d:/data_analysis/zls_project/0SPD_dual-sg_pertubseq_project/2025-Nature_genetics_respone_prepare/2025-6-6-substype_analysis_In_Ex_neuron/")
# Load required packages
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(rmarkdown) #to get rmd file 

# Read and prepare the data
# for Ex_neuron
data <- read.csv("Ex_neuron.MDD.riskGene.DEGs.in.PseudoBulk.csv") %>%
  mutate(Per_genes = as.character(Per_genes),
         layer = as.character(layer),
         gene = as.character(gene))
         
# for In_neuron
#data <- read.csv("In_neuron.MDD.riskGene.DEGs.in.PseudoBulk.csv") %>%
#  mutate(Per_genes = as.character(Per_genes),
#         layer = as.character(layer),
#         gene = as.character(gene))

# Filter for significant DEGs (FDR < 0.05)
sig_data <- data %>% 
  filter(FDR < 0.05) %>%
  select(layer, Per_genes, gene, logFC)

# Identify genes unique to each Per_gene-layer combination
specific_genes <- sig_data %>%
  group_by(gene) %>%
  filter(n() == 1) %>%  # Genes that appear only once (unique to one Per_gene-layer combo)
  ungroup()

# Create a combined identifier for rows (Per_gene + layer)
specific_genes <- specific_genes %>%
  unite("Per_layer", Per_genes, layer, sep = "|", remove = FALSE)

# Get all unique specific genes and Per_gene-layer combinations
all_genes <- unique(specific_genes$gene)
all_per_layers <- unique(specific_genes$Per_layer)

# Create a matrix for the heatmap
heatmap_matrix <- matrix(NA, 
                         nrow = length(all_per_layers),
                         ncol = length(all_genes),
                         dimnames = list(all_per_layers, all_genes))

# Fill the matrix with logFC values
for (i in 1:nrow(specific_genes)) {
  row_idx <- which(rownames(heatmap_matrix) == specific_genes$Per_layer[i])
  col_idx <- which(colnames(heatmap_matrix) == specific_genes$gene[i])
  heatmap_matrix[row_idx, col_idx] <- specific_genes$logFC[i]
}

# Prepare annotation data
layer_info <- str_split_fixed(rownames(heatmap_matrix), "\\|", 2)[,2]
per_gene_info <- str_split_fixed(rownames(heatmap_matrix), "\\|", 2)[,1]
library(RColorBrewer)
# Create annotations
n_layers <- length(unique(layer_info))
n_per_genes <- length(unique(per_gene_info))

# bulid color for Layer
if(n_layers <= 8) {
  layer_colors <- setNames(brewer.pal(n_layers, "Set2"), unique(layer_info))
} else if(n_layers <= 12) {
  layer_colors <- setNames(brewer.pal(n_layers, "Set3"), unique(layer_info))
} else {
  layer_colors <- setNames(colorRampPalette(brewer.pal(8, "Set2"))(n_layers), unique(layer_info))
}

# build color for Per_gene
if(n_per_genes <= 8) {
  per_gene_colors <- setNames(brewer.pal(n_per_genes, "Dark2"), unique(per_gene_info))
} else if(n_per_genes <= 12) {
  per_gene_colors <- setNames(brewer.pal(n_per_genes, "Paired"), unique(per_gene_info))
} else {
  per_gene_colors <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(n_per_genes), unique(per_gene_info))
}

# build row for annotation
row_ha <- rowAnnotation(
  Layer = layer_info,
  Per_gene = per_gene_info,
  col = list(
    Layer = layer_colors,
    Per_gene = per_gene_colors
  ),
  show_legend = c(TRUE, TRUE),
  show_annotation_name = TRUE,
  annotation_legend_param = list(
    Layer = list(title = "Layer"),
    Per_gene = list(title = "Per_gene")
  )
)


# Create the heatmap
ht <- Heatmap(heatmap_matrix,
              name = "logFC",
              col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")),
              na_col = "gray90",
              show_row_names = FALSE,
              show_column_names = TRUE,
              column_names_gp = gpar(fontsize = 4),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              right_annotation = row_ha,
              heatmap_legend_param = list(title = "logFC"),
              column_title = "Genes specifically perturbed by each Per_gene in each layer")
# Draw the heatmap
draw(ht, heatmap_legend_side = "right", 
     column_title = "Perturbation Screening Results", 
     column_title_gp = gpar(fontsize = 14))

# Save the results
#for Ex_neuron
write.csv(specific_genes, "Per_gene_layer_specific_DEGs.Ex_neuron.csv", row.names = FALSE)

pdf("Per_gene_specific_DEGs_heatmap_with_layers.Ex_neuron.pdf", width = 15.68, height = 7.75)
draw(ht, heatmap_legend_side = "right", 
     column_title = "Perturbation Screening Results", 
     column_title_gp = gpar(fontsize = 14))
dev.off()

#render("Perturbation_Results.Ex_neuron.Rmd", output_file = "Perturbation_Results.Ex_neuron.html")



#for In_neuron

#write.csv(specific_genes, "Per_gene_layer_specific_DEGs.In_neuron.csv", row.names = FALSE)

#pdf("Per_gene_specific_DEGs_heatmap_with_layers.In_neuron.pdf", width = 15.68, height = 7.75)
#draw(ht, heatmap_legend_side = "right", 
#     column_title = "Perturbation Screening Results", 
#     column_title_gp = gpar(fontsize = 14))
#dev.off()

#render("Perturbation_Results.In_neuron.Rmd", output_file = "Perturbation_Results.In_neuron.html")

