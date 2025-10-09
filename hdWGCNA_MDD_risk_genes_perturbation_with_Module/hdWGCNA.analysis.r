# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)
library(devtools)
library(igraph)

# co-expression network analysis packages:
library(WGCNA)
# devtools::install_github('smorabit/hdWGCNA', ref='dev')
library(hdWGCNA)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)
seurat_obj <- readRDS('../../analysis20250219/sPD.deBatch.merged_final2.1moi.fix.rds') 
#the input file you can download from OMIX database of NGDC (see OMIX010695: https://ngdc.cncb.ac.cn/omix/preview/msLVjcbG)
#head(seurat_obj@meta.data)
#                                                 cellname orig.ident nCount_RNA
#ACGATGTAGTGGAGTC-1_1_1_1_1_1 ACGATGTAGTGGAGTC-1_1_1_1_1_1      SPD_2       6961
#ACGATGTAGTGGCACA-1_1_1_1_1_1 ACGATGTAGTGGCACA-1_1_1_1_1_1      SPD_2       5317
#ACGATGTCAATAGAGT-1_1_1_1_1_1 ACGATGTCAATAGAGT-1_1_1_1_1_1      SPD_2      27059
#ACGATGTCACCGCTAG-1_1_1_1_1_1 ACGATGTCACCGCTAG-1_1_1_1_1_1      SPD_2      16369
#ACGATGTCATGGTTGT-1_1_1_1_1_1 ACGATGTCATGGTTGT-1_1_1_1_1_1      SPD_2      19177
#ACGATGTGTCCAAGTT-1_1_1_1_1_1 ACGATGTGTCCAAGTT-1_1_1_1_1_1      SPD_2      11683
#                             nFeature_RNA  percent.mt integrated_snn_res.0.2
#ACGATGTAGTGGAGTC-1_1_1_1_1_1         3288 0.000000000                      0
#ACGATGTAGTGGCACA-1_1_1_1_1_1         2579 0.000000000                     16
#ACGATGTCAATAGAGT-1_1_1_1_1_1         6882 0.003695628                     17
#ACGATGTCACCGCTAG-1_1_1_1_1_1         4377 0.000000000                     13
#ACGATGTCATGGTTGT-1_1_1_1_1_1         5586 0.015643740                     10
#ACGATGTGTCCAAGTT-1_1_1_1_1_1         4078 0.000000000                      1
#                             seurat_clusters integrated_snn_res.0.3 cell_anno
#ACGATGTAGTGGAGTC-1_1_1_1_1_1               0                      0        t1
#ACGATGTAGTGGCACA-1_1_1_1_1_1              16                     16       t17
#ACGATGTCAATAGAGT-1_1_1_1_1_1              18                     18       t19
#ACGATGTCACCGCTAG-1_1_1_1_1_1              14                     14       t15
#ACGATGTCATGGTTGT-1_1_1_1_1_1              10                     10       t11
#ACGATGTGTCCAAGTT-1_1_1_1_1_1               1                      1        t2
#                              celltype  celltype2            cellBarcode moi
#ACGATGTAGTGGAGTC-1_1_1_1_1_1 In_neuron     neuron SPD_2_ACGATGTAGTGGAGTC   1
#ACGATGTAGTGGCACA-1_1_1_1_1_1 In_neuron     neuron SPD_2_ACGATGTAGTGGCACA   1
#ACGATGTCAATAGAGT-1_1_1_1_1_1 Microglia non_neuron SPD_2_ACGATGTCAATAGAGT   1
#ACGATGTCACCGCTAG-1_1_1_1_1_1 Ex_neuron     neuron SPD_2_ACGATGTCACCGCTAG   1
#ACGATGTCATGGTTGT-1_1_1_1_1_1 Ex_neuron     neuron SPD_2_ACGATGTCATGGTTGT   1
#ACGATGTGTCCAAGTT-1_1_1_1_1_1 Ex_neuron     neuron SPD_2_ACGATGTGTCCAAGTT   1
#                                      id          sg2_rev_seq
#ACGATGTAGTGGAGTC-1_1_1_1_1_1   zls_Erbb4 GCATCGAGCACAACCGGGAC
#ACGATGTAGTGGCACA-1_1_1_1_1_1  zls_Ccser1 ATGTCAGCTGCATGCCCAAC
#ACGATGTCAATAGAGT-1_1_1_1_1_1 zls_Dennd1a CTGGGAATTGCCTCTGCACC
#ACGATGTCACCGCTAG-1_1_1_1_1_1   zls_Cnnm2 CATTGTTCACTCGCTGCACA
#ACGATGTCATGGTTGT-1_1_1_1_1_1   zls_Erbb4 GCATCGAGCACAACCGGGAC
#ACGATGTGTCCAAGTT-1_1_1_1_1_1    zls_Snrk ACGACAGCGAGACACTGACC
#                                      sg1_rev_seq sg2_read_count sg2_umi_count
#ACGATGTAGTGGAGTC-1_1_1_1_1_1 CGCTGCCATCAGAAGGCTCC           2984            35
#ACGATGTAGTGGCACA-1_1_1_1_1_1 CTGGTGTGCATAGCTCCTCC            659            13
#ACGATGTCAATAGAGT-1_1_1_1_1_1 GTCCTAGGATAAGCCACTTC           9115           186
#ACGATGTCACCGCTAG-1_1_1_1_1_1 ACGCTCGCTGGGTGTGTGCC           1554            21
#ACGATGTCATGGTTGT-1_1_1_1_1_1 CGCTGCCATCAGAAGGCTCC           4673            64
#ACGATGTGTCCAAGTT-1_1_1_1_1_1 AGTCTTCATTAAGACCCTCC            768            15
#                             sg2_coverage sg2_coverage_threshold sg2_proportion
#ACGATGTAGTGGAGTC-1_1_1_1_1_1     85.25714                   TRUE      0.4375000
#ACGATGTAGTGGCACA-1_1_1_1_1_1     50.69231                   TRUE      0.2826087
#ACGATGTCAATAGAGT-1_1_1_1_1_1     49.00538                   TRUE      0.3496241
#ACGATGTCACCGCTAG-1_1_1_1_1_1     74.00000                   TRUE      0.4200000
#ACGATGTCATGGTTGT-1_1_1_1_1_1     73.01562                   TRUE      0.2126246
#ACGATGTGTCCAAGTT-1_1_1_1_1_1     51.20000                   TRUE      0.2142857
#                             sg2_moi sg1_read_count sg1_umi_count sg1_coverage
#ACGATGTAGTGGAGTC-1_1_1_1_1_1       1           2282            38     60.05263
#ACGATGTAGTGGCACA-1_1_1_1_1_1       1            141             4     35.25000
#ACGATGTCAATAGAGT-1_1_1_1_1_1       1            208             3     69.33333
#ACGATGTCACCGCTAG-1_1_1_1_1_1       1           1846            41     45.02439
#ACGATGTCATGGTTGT-1_1_1_1_1_1       1           1283            18     71.27778
#ACGATGTGTCCAAGTT-1_1_1_1_1_1       1           2470            52     47.50000
#                             sg1_coverage_threshold sg1_proportion
#ACGATGTAGTGGAGTC-1_1_1_1_1_1                   TRUE    0.316666667
#ACGATGTAGTGGCACA-1_1_1_1_1_1                   TRUE    0.100000000
#ACGATGTCAATAGAGT-1_1_1_1_1_1                   TRUE    0.008902077
#ACGATGTCACCGCTAG-1_1_1_1_1_1                   TRUE    0.621212121
#ACGATGTCATGGTTGT-1_1_1_1_1_1                   TRUE    0.057877814
#ACGATGTGTCCAAGTT-1_1_1_1_1_1                   TRUE    0.481481481
#write.table(seurat_obj@meta.data,"sPD.deBatch.merged_final2.1moi.fix.metadata.txt",quote=F,row.names=F) 
#Define cells and samples to ensure that the analysis can be conducted strictly according to the hdWGCNA protocol; **each MDD risk gene is considered as a sample**.
seurat_obj@meta.data$cell_type<-seurat_obj@meta.data$celltype
seurat_obj@meta.data$cell_type2<-seurat_obj@meta.data$celltype2
seurat_obj@meta.data$Sample<-seurat_obj@meta.data$id
seurat_obj@meta.data[grepl("WBL",seurat_obj@meta.data$id),]$Sample<-"control"
seurat_obj@meta.data[grepl("SPD87",seurat_obj@meta.data$id),]$Sample<-"control"
seurat_obj@meta.data[grepl("zls_ctrl",seurat_obj@meta.data$id),]$Sample<-"control"


#table(seurat_obj@meta.data$Sample)
#> table(seurat_obj@meta.data$Sample)

#     control   zls_Acvr1b    zls_Agbl4    zls_Ano10    zls_Ap3b1  zls_Arfgef2
#        2277           71           79            5         1802           41
#   zls_Ascc3    zls_Astn1    zls_Astn2   zls_Atp2a2 zls_B4galnt4    zls_Baz2b
#         126           30         1531           25          812          202
#   zls_Bltp1   zls_Ccser1     zls_Cdh9    zls_Chmp3    zls_Cnnm2    zls_Cnot4
#         129         3981          331          907         1544           53
#   zls_Cntln   zls_Cstpp1  zls_Cttnbp2    zls_Dagla  zls_Dennd1a  zls_Dennd1b
#         252         1144           83          171         1729           48
#  zls_Elavl2   zls_Elavl4    zls_Erbb4    zls_Esrrg     zls_Etv1     zls_Etv6
#          48           36         1980           18          173            9
#   zls_Exoc4  zls_Fam120a   zls_Fbxl17     zls_Fhit    zls_Foxp2   zls_Gabra1
#         612          286           22          232          138           71
#  zls_Galnt2   zls_Gigyf2     zls_Grm5      zls_Htt    zls_Kcnq5  zls_Kirrel3
#         669          216           45           41          234          136
#    zls_Lhpp    zls_Lrp1b    zls_Lsamp    zls_Ltbp3    zls_Luzp2    zls_Maml3
#         337           31          114          181          580           19
#    zls_Nbas   zls_Nkain2    zls_Nlgn1    zls_Nrxn1    zls_Pacrg     zls_Pclo
#          43         2350          810          232          229          167
#   zls_Plcl2    zls_Qser1     zls_Rbks     zls_Rere    zls_Rev3l    zls_Samd5
#         841           51          608          446           60           61
#    zls_Sdk1     zls_Sgcz    zls_Sgip1   zls_Shank2 zls_Slc25a12  zls_Slc30a9
#         107           44            4           23          103          237
#    zls_Snrk   zls_Sorcs3     zls_Sox5    zls_Sppl3    zls_Suds3     zls_Tank
#        1979           70           52          360         1672          184
#    zls_Tcf4    zls_Tenm2      zls_Tox     zls_Usp3   zls_Zfp638  zls_Zfp804a
#          19          144          715          159           22          621
#  zls_Zmynd8
#          23
###
p <- DimPlot(seurat_obj, group.by='cell_type', label=TRUE)+
  ggtitle('Neuron cells of pool sgRNA screening of MDD risk genes') + NoLegend()
pdf("f1.pdf")
p
dev.off()
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "tutorial" # the name of the hdWGCNA experiment
)

   #length(seurat_obj@misc$tutorial$wgcna_genes)
   #table(seurat_obj$cell_type)
   #table(seurat_obj$Sample)
   #seurat_obj

# construct metacells  in each group
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("cell_type", "Sample"), # specify the columns in seurat_obj@meta.data to group by
  k = 25, # nearest-neighbors parameter
  reduction = "umap", 
  slot='counts',
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'cell_type' # set the Idents of the metacell seurat object
)

seurat_obj <- NormalizeMetacells(seurat_obj)

  #  get the metacell object from the hdWGCNA experiment
  #  metacell_obj <- GetMetacellObject(seurat_obj)

  #seurat_obj
  # table(seurat_obj$cell_type,seurat_obj$Sample)
  # table(metacell_obj$cell_type,metacell_obj$Sample)
 
 
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = c("In_neuron","Ex_neuron"), # the name of the group of interest in the group.by column
  group.by='cell_type', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  use_metacells = TRUE, # use the metacells (TRUE) or the full expression matrix (FALSE)
  slot = 'data' # using normalized data
)

# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
   powers = c(seq(1, 10, by = 1), seq(12, 30, by = 2))) 
  # networkType = 'signed' 

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj,
                            point_size = 5,
                            text_size = 3)
                            # selected_power = NULL
# assemble with patchwork
pdf('f2.pdf')
wrap_plots(plot_list, ncol=2)
dev.off()


power_table <- GetPowerTable(seurat_obj)
write.table(power_table,"hdWGCNA.power_table.txt",quote=F,row.names=F)

#after check f2.pdf£¬we get soft_power number is 6

seurat_obj <- ConstructNetwork(
  seurat_obj,
  soft_power=6,
  setDatExpr=FALSE,
  corType = "pearson",
  networkType = "signed",
  TOMType = "signed",
  detectCutHeight = 0.995,
  minModuleSize = 50,
  mergeCutHeight = 0.2,
  tom_outdir = "TOM", # output directory
  tom_name = 'Neuron' # name of the topoligical overlap matrix written to disk
)

pdf('f3.pdf')  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>produce figure 2a
PlotDendrogram(seurat_obj, main='Neuron hdWGCNA Dendrogram')
dev.off()

# need to run ScaleData first or else harmony throws an error:
seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  scale.model.use = "linear", #  choices are "linear", "poisson", or "negbinom"
  assay = NULL, # DefaultAssay(seurat_obj)
  pc_dim = 1,
  group.by.vars="Sample" # harmonize
)

########
cur_traits <- c('sg2_coverage','sg1_coverage','integrated_snn_res.0.3',
                'nCount_RNA', 'nFeature_RNA', 'percent.mt')
seurat_obj <- ModuleTraitCorrelation(
  seurat_obj,
  traits = cur_traits,
  features = "hMEs", # Valid choices are hMEs, MEs, or scores
  cor_method = "pearson", # Valid choices are pearson, spearman, kendall.
  group.by='cell_type'
)


# get the mt-correlation results
mt_cor <- GetModuleTraitCorrelation(seurat_obj)

names(mt_cor)
# "cor"  "pval" "fdr"


pdf('f4.pdf')
PlotModuleTraitCorrelation(
  seurat_obj,
  label = 'fdr', # add p-val label in each cell of the heatmap
  label_symbol = 'stars', # labels as 'stars' or as 'numeric'
  text_size = 2,
  text_digits = 2,
  text_color = 'white',
  high_color = '#fc9272',
  mid_color = '#ffffbf',
  low_color = '#9ecae1',
  plot_max = 0.2,
  combine=T 
)

PlotModuleTraitCorrelation(
  seurat_obj,
  label = 'fdr', # add p-val label in each cell of the heatmap
  label_symbol = 'stars', # labels as 'stars' or as 'numeric'
  text_size = 2,
  text_digits = 2,
  text_color = 'white',
  high_color = '#fc9272',
  mid_color = '#ffffbf',
  low_color = '#9ecae1',
  plot_max = 0.2,
  combine=F 
)

dev.off()


# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'cell_type',
  corFnc = "bicor", # to obtain Pearson correlation
  corOptions = "use='p'", # to obtain Pearson correlation
  harmonized = TRUE,
  assay = NULL,
  slot = "data", # default to normalized 'data' slot
  group_name = c('In_neuron','Ex_neuron') # select only neuron cells
)
# rename the modules 
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "Neuron-M"
)

modules <- GetModules(seurat_obj)
write.table(modules,"hdWGCNA.modules.txt",quote=F)

# plot genes ranked by kME for each module
p <- PlotKMEs(seurat_obj,
              ncol=5,
              n_hubs = 10, # number of hub genes to display
              text_size = 2,
              plot_widths = c(3, 2) # the relative width between the kME rank plot and the hub gene text
              )
pdf('f5.pdf') #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>produce figure 2b
p
dev.off()

# get hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)
write.table(hub_df,"hdWGCNA.hub_df.txt",quote=F,row.names=F)
#
#saveRDS(seurat_obj, file='hdWGCNA_object.tmp.rds')

library(UCell)
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell'
)

##################################################
# make a featureplot of hMEs for each module
plot_list <-ModuleFeaturePlot(
  seurat_obj,
  reduction = "umap",
  features = "hMEs",
  order_points = TRUE, # order so the points with highest hMEs are on top
  restrict_range = TRUE,
  point_size = 0.5,
  alpha = 1,
  label_legend = FALSE,
  raster_dpi = 500,
  raster_scale = 1,
  plot_ratio = 1,
  title = TRUE
)
# stitch together with patchwork
pdf("moduleFeaturePlot.hMEs.pdf")
wrap_plots(plot_list, ncol=5)
dev.off()


# make a featureplot of hub scores for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='scores', # plot the hub gene scores
  order='shuffle', # order so cells are shuffled
  ucell = TRUE # depending on Seurat vs UCell for gene scoring
)

# stitch together with patchwork
pdf("moduleFeaturePlot.scores.pdf")
wrap_plots(plot_list, ncol=5)
dev.off()
##################################################

# plot module correlagram
pdf("moduleFeaturePlot.correlation.pdf") #
# plot module correlagram
ModuleCorrelogram(seurat_obj,
                  exclude_grey = TRUE, 
                  features = "hMEs" # What to plot? Can select hMEs, MEs, scores, or average
                  )
# plot module correlagram
ModuleCorrelogram(seurat_obj,
                  exclude_grey = TRUE, 
                  features = "MEs" # What to plot? Can select hMEs, MEs, scores, or average
                  )
# plot module correlagram£¬#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>produce Extend Data Fig.6c
ModuleCorrelogram(seurat_obj,
                  exclude_grey = TRUE, 
                  features = "scores" # What to plot? Can select hMEs, MEs, scores, or average
                  )
dev.off()


# get hMEs from seurat object
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
mods <- colnames(MEs); 
mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

# plot with Seurat's DotPlot function
p <- DotPlot(seurat_obj, features=mods, group.by = 'cell_type')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

#Image
pdf("bulbPlot.pdf")
p
dev.off()

plot_list <- lapply(mods, function(x) {
  print(x)
  p <- VlnPlot(
    seurat_obj,
    features = x,
    group.by = 'cell_type',
    pt.size = 0 # don't show actual data points
  )

  # add box-and-whisker plots on top:
  p <- p + geom_boxplot(width=.25, fill='white')

  # change axis labels and remove legend:
  p <- p + xlab('') + ylab('hME') + NoLegend()

  p
})

pdf("hME.all.pdf",width=20,height=10)
wrap_plots(plot_list, ncol = 5)
dev.off()


library(igraph)
# Visualizes the top hub genes for selected modules as a circular network plot
ModuleNetworkPlot(
  seurat_obj,
  mods = "all", # all modules are plotted.
  outdir = "ModuleNetworks", # The directory where the plots will be stored.
  plot_size = c(6, 6),
  label_center = FALSE,
  edge.alpha = 0.25,
  vertex.label.cex = 1, 
  vertex.size = 6 
)
# hubgene network
options(future.globals.maxSize = 10 * 1024^3)  # 10 GB

pdf("hubGeneNetwork.pdf")
HubGeneNetworkPlot(
  seurat_obj,
  mods = "all", # all modules are plotted.
  n_hubs = 3,
  n_other=6,
  edge_prop = 0.75,
)
dev.off()

seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 10, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1 # min distance between points in UMAP space
)
# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)

pdf("modules.umap.pdf") #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>produce Extend Data Fig.6b
# get the hub gene UMAP table from the seurat object
# plot with ggplot
ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
   color=umap_df$color, # color each point by WGCNA module
   size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()

#
ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, 
    size=umap_df$kME*2  
  ) +
  scale_x_continuous(name="UMAP1") +
  scale_y_continuous(name="UMAP2") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(color = "black", size = 12), 
    axis.title.y = element_text(color = "black", size = 12), 
    axis.text.x = element_text(color = "black", size = 10),  
    axis.text.y = element_text(color = "black", size = 10)   
  ) 
dev.off()

dev.off()

pdf("modules.umap2.pdf")
ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=2 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE

)
dev.off()

pdf("RadarPlot.pdf")
ModuleRadarPlot(
  seurat_obj,
  group.by ='seurat_clusters',
  barcodes = seurat_obj@meta.data %>% subset(celltype2 =='neuron') %>% rownames(),
  axis.label.size=4,
  grid.label.size=4
)

saveRDS(seurat_obj,"hdWGCNA_object.final.rds")
#
#seurat_obj<-readRDS("hdWGCNA_object.final.rds")
#
#next part enrichment analysis for each module
# define the enrichr databases to test
dbs <- c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023')

# perform enrichment tests
seurat_obj <- RunEnrichr(
  seurat_obj,
  dbs=dbs,
  max_genes = 100 # use max_genes = Inf to choose all genes
)

# retrieve the output table
enrich_df <- GetEnrichrTable(seurat_obj)
#saveRDS(enrich_df,"enrich_df.rds")
#enrich_df <-readRDS("enrich_df.rds")
##plot by home code
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(purrr)

module_enrich_df<-enrich_df
modules <- unique(module_enrich_df$module)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>produce Extend Data Fig.6d-6m
output_dir <- "Enrichr_plot3"
if (!dir.exists(output_dir)) dir.create(output_dir)

top10_by_module <- module_enrich_df %>%
  mutate(Term_simple = Term ) %>%  
  arrange(P.value) %>%
  group_by(module) %>%
  slice_min(P.value, n = 10) %>%
  ungroup() %>%
  mutate(
    log10P = -log10(P.value),
    log10CombinedScore = log10(Combined.Score)
  ) %>%
  split(.$module)

plot_top10 <- function(data, mod_name) {
  ggplot(data, aes(x = fct_reorder(Term_simple, log10P), 
                   y = log10P, 
                   fill = log10CombinedScore)) +
    geom_bar(stat = "identity", width = 0.8) +
    coord_flip() +
    scale_fill_gradient(low = "lightblue", high = "darkblue") +
    labs(title = paste("Module:", mod_name),
         x = NULL, y = "-log10(P-value)") +
    theme_minimal()
}

walk2(top10_by_module, names(top10_by_module), ~ {
  ggsave(
    filename = file.path(output_dir, paste0("module_NEW.", gsub("\\W", "_", .y), ".pdf")),
    plot = plot_top10(.x, .y),
    width = 8, height = 6
  )
})


