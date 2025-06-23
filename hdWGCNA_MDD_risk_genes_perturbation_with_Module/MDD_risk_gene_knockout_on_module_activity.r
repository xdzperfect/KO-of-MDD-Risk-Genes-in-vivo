#2025-6-13
#figure for MDD risk genes correlative to hdWGCNA module
rm(list=ls())
# single-cell analysis package
library(Seurat)
library(circlize)
library(magick)
library(dendsort)
library(dplyr)
library(ComplexHeatmap)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("WGCNA")
#install.packages("hdWGCNA-dev.zip")
# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 8)

setwd("d:/data_analysis/zls_project/0SPD_dual-sg_pertubseq_project/2025-Nature_genetics_respone_prepare/2025-5-12-hdWGCNA-MDD/")
seurat_obj<-readRDS("hdWGCNA_object.final2.rds") #the file you can download from the National Genomics Data Center (NGDC) with a project ID of PRJCA031709, in OMIX(OMIX010695)
seurat_obj@misc$tutorial$wgcna_net$TOMFiles<-"D:/data_analysis/zls_project/0SPD_dual-sg_pertubseq_project/2025-Nature_genetics_respone_prepare/2025-5-12-hdWGCNA-MDD/TOM/Neuron_TOM.rda"
#the file you can download from the National Genomics Data Center (NGDC) with a project ID of PRJCA031709, in OMIX(OMIX010695)
#head(seurat_obj@meta.data)
module_df <- seurat_obj@meta.data %>% 
  select(starts_with("Neuron-M"))  
module_colors <- colnames(module_df) 
module_colors_sorted <- module_colors[order(as.numeric(gsub("Neuron-M", "", module_colors)))]
#str(seurat_obj@meta.data)
#table(seurat_obj@meta.data$id)
seurat_obj@meta.data$perturbation<-gsub("zls_","",seurat_obj@meta.data$id)
seurat_obj@meta.data[grepl("WBL|SPD|ctrl",seurat_obj@meta.data$id),]$perturbation<-"control"

ko_module_matrix <- aggregate(
  module_df,
  by = list(perturbation = seurat_obj@meta.data$perturbation),  # 替换为实际的KO列名
  FUN = mean
) %>%
  tibble::column_to_rownames("perturbation")

ko_module_matrix <- ko_module_matrix[, module_colors_sorted]

control_means <- ko_module_matrix["control", ]
ko_matrix <- ko_module_matrix[rownames(ko_module_matrix) != "control", ]
ko_matrix <- as.matrix(ko_matrix)
control_means <- as.numeric(control_means)
ko_diff_matrix <- sweep(ko_matrix, 2, control_means, "-")
library(pheatmap)
p1<-pheatmap(
  ko_diff_matrix,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-3, 3, length.out = 100),  
  main = "KO Effects Relative to Control",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  fontsize_row = 8,
  fontsize_col = 8
)

pdf("KO_MDD_risk_Gene_to_Module1.pdf",width =7.20 ,height = 9)
p1
dev.off()


