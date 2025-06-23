#(base) xuzhzh@xdzpefect:/mnt/disk12_part2/zls_project/publicData_SPD/seuratV5_query_test/mdd_query_raw$ cat 2025-6-11.r
#ref data (MDD patients)
library(Seurat)
library(ggplot2)

pancreas.ref<-readRDS("../deBatch.merged.resolution0.1.rds")
#the data you can download from NGDC database-project ID: PRJCA031709
meta<-pancreas.ref@meta.data
ref<-read.table("../GSE213982_combined_counts_matrix_cells_columns.txt",header=F)
ref<-ref[grepl("F",ref$V1),]
colnames(ref)<-c("cellID","celltype","celltype_sub")
meta$cellID<-row.names(meta)
meta<-merge(meta,ref,by="cellID")
pancreas.ref@meta.data$cellID<-row.names(pancreas.ref@meta.data)
meta2<-meta[match(pancreas.ref@meta.data$cellID,meta$cellID),]
pancreas.ref@meta.data<-meta2
row.names(pancreas.ref@meta.data)<-pancreas.ref@meta.data$cellID
#query data (pool screening of KO MDD risk Genes)
merged<-readRDS("/mnt/disk12_part2/zls_project/sPD_project/snRNA-seq/single2024-9-4/resolution0.3/analysis20250219/sPD.deBatch.merged_final2.1moi.fix.rds")
#the data you can download from NGDC database-project ID: PRJCA031709
#
ko_data<-CreateSeuratObject(counts=merged[["RNA"]]$counts,project='ko_query',min.cells=0,min.features=0)
ko_data
ko_data@meta.data$celltype<-merged@meta.data$celltype
pancreas.query <- NormalizeData(ko_data)
pancreas.query@meta.data$celltype<-gsub("_neuron","N",pancreas.query@meta.data$celltype)
pancreas.query@meta.data$celltype<-gsub("Astrocyte","Ast",pancreas.query@meta.data$celltype)
pancreas.query@meta.data$celltype<-gsub("Microglia","Mic",pancreas.query@meta.data$celltype)

head(pancreas.ref@meta.data)
#find Anchors and MapQuery

#get ref data UMAP value
pancreas.ref <- RunUMAP(pancreas.ref, dims = 1:30, reduction = "pca", return.model = TRUE)
pancreas.anchors <- FindTransferAnchors(reference = pancreas.ref, query = pancreas.query, dims = 1:30, reference.reduction = "pca", k.anchor=30)
#Projecting cell embeddings
#Finding neighborhoods
#Finding anchors
#        Found 3726 anchors
pancreas.query <- MapQuery(anchorset = pancreas.anchors, reference = pancreas.ref, query = pancreas.query, refdata = list(celltype = "celltype"), reference.reduction = "pca", reduction.model ="umap")
predictions <- TransferData(anchorset = pancreas.anchors, refdata = pancreas.ref$celltype, dims = 1:30)
pancreas.query <- AddMetaData(pancreas.query, metadata = predictions)
pancreas.query$prediction.match <- pancreas.query$predicted.id == pancreas.query$celltype

table(pancreas.query$prediction.match)
#FALSE  TRUE
#18154 17883

table(pancreas.query$predicted.id)
#  Ast   ExN   Mic
#  923 34583   531

pdf("vlnPlot.tmp3b.pdf")

VlnPlot(pancreas.query, c("Aqp4", "Gja1", "Rbfox3", "Grin1", "C1qa", "C1qb"), group.by = "predicted.id")
dev.off()
table(pancreas.query@meta.data$celltype)
#              Ast Endothelial_cells               ExN               InN
#             1130               134             17572             16645
#              Mic              OPCs   Oligodendrocyte
#              328                83               145


p1<-DimPlot(pancreas.ref, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 3,repel = TRUE) + NoLegend() + ggtitle("Reference annotations") +coord_cartesian(xlim = c(-15, 15),ylim=c(-20,20))
p2<-DimPlot(pancreas.query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels") +coord_cartesian(xlim = c(-15, 15),ylim=c(-20,20))

pdf("umap.tmp2b.pdf",width=20,height=10)
#pdf("umap.tmp2b.2025-6-22.pdf",width=20,height=10)
p1+p2
dev.off()

pancreas.query@meta.data$id<-gsub("zls_","",merged@meta.data$id)
pancreas.query@meta.data$orig.ident<-gsub("zls_","",merged@meta.data$orig.ident)
pancreas.query@meta.data$cellBarcode<-merged@meta.data$cellBarcode
pancreas.query@meta.data$moi<-merged@meta.data$moi
pancreas.query$treatment<-"treatment"
pancreas.query@meta.data[grepl("WBL|SPD|ctrl",pancreas.query@meta.data$id),]$treatment<-"control"
table(pancreas.query@meta.data$treatment)
# control treatment
     2277     33760
saveRDS(pancreas.query,"sPD.1moi.query2ref_MDD.rds")

#get only the cells of sPD project to MDD
Idents(pancreas.query)<-pancreas.query@meta.data[,"prediction.match"]
sub_match<-subset(pancreas.query,idents="TRUE")
saveRDS(sub_match,"sPD.1moi.query2ref_MDD.matched.rds")
#
Idents(sub_match)<-sub_match@meta.data[,"predicted.celltype"]
Ex_neuron<-subset(sub_match,idents="ExN")
saveRDS(Ex_neuron,"sPD.1moi.query2ref_MDD.matched.Ex_neuron.rds")

