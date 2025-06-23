#MDD patients reanalysis
cd /mnt/disk12_part2/zls_project/publicData_SPD/resolution0.1
#>>>>>>>>>>>>>>>>>>>>>>cell annotaion<<<<<<<<<<<<<<<<<<<<<<<<

(R4.3.1) xuzhzh@xdzpefect:/mnt/disk12_part2/zls_project/publicData_SPD/resolution0.1$ cat resolution0.1.r
library("Seurat")
merged<-readRDS("../deBatch.merged.resolution0.1.rpca.rds")
#merged.markers <- FindAllMarkers(merged, min.pct = 0.25, logfc.threshold = 0.25)
#write.csv(merged.markers,"sPD.deBatch.merged.markers.csv")

new.cluster.ids <- c("t1","t2","t3","t4","t5","t6","t7","t8","t9","t10","t11","t12","t13","t14","t15","t16","t17","t18")
names(new.cluster.ids) <- levels(merged)

merged <- RenameIdents(merged, new.cluster.ids)
pdf("dimPlot.debatch.addCellNames.pdf")
DimPlot(merged, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(merged, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

merged@meta.data$cell_anno <- Idents(merged)
write.csv(merged@meta.data, file="sPD.deBatch.metadata.merged.csv")
saveRDS(merged,file="sPD.deBatch.merged_final.rds")
#merged<-readRDS("sPD.deBatch.merged_final.rds")

library(celldex)
library(SingleR)
##mouse
#load("/mnt/disk11_8T/project_maq_aging/snRNA-seq/resolution0.1_debatch/mouse.sc.RData")
#ref<-mouse.sc
#human 
load("/mnt/disk10_8T/zls_project/public_data/PD/GSE202210/human.sc.RDdata")
ref<-human.sc
#
refdata<-ref
#
testdata <- GetAssayData(merged,slot = "data")
clusters <- merged@meta.data$seurat_clusters
cellpred <- SingleR( test=testdata,
                     ref = refdata,
                     labels = refdata$label.main,
                     method = "cluster",
                     clusters = clusters,
                     assay.type.test = "logcounts",
                     assay.type.ref = "logcounts")
celltype = data.frame(ClusterID = rownames(cellpred),
                      celltype = cellpred$labels,
                      stringsAsFactors = F)
#celltype

write.csv(celltype,"celltype_anno_SingleR.resolution0.1.csv")
save.image("2024-4-1.singleR.RData")
pdf("plotScoreHeatmap.resolution0.1.pdf")
p<- plotScoreHeatmap(cellpred,clusters = rownames(cellpred),
                     order.by = "cluster")
p
dev.off()
#


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>add cell type and get cluster1 genes' expression data part<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
cat >2024-9-30.get_female_and_allPatient.r
library(Seurat)
#for all (Ex_N and IN_N)
merged<-readRDS("sPD.deBatch.merged_final.rds")
merged@meta.data$sampleType<-"WT"
merged@meta.data[merged@meta.data$orig.ident=="F1",]$sampleType<-"MDD"
merged@meta.data[merged@meta.data$orig.ident=="F11",]$sampleType<-"MDD"
merged@meta.data[merged@meta.data$orig.ident=="F12",]$sampleType<-"MDD"
merged@meta.data[merged@meta.data$orig.ident=="F14",]$sampleType<-"MDD"
merged@meta.data[merged@meta.data$orig.ident=="F15",]$sampleType<-"MDD"
merged@meta.data[merged@meta.data$orig.ident=="F16",]$sampleType<-"MDD"
merged@meta.data[merged@meta.data$orig.ident=="F17",]$sampleType<-"MDD"
merged@meta.data[merged@meta.data$orig.ident=="F18",]$sampleType<-"MDD"
merged@meta.data[merged@meta.data$orig.ident=="F19",]$sampleType<-"MDD"
merged@meta.data[merged@meta.data$orig.ident=="F2",]$sampleType<-"MDD"
merged@meta.data[merged@meta.data$orig.ident=="F20",]$sampleType<-"MDD"
merged@meta.data[merged@meta.data$orig.ident=="F25",]$sampleType<-"MDD"
merged@meta.data[merged@meta.data$orig.ident=="F27",]$sampleType<-"MDD"
merged@meta.data[merged@meta.data$orig.ident=="F28",]$sampleType<-"MDD"
merged@meta.data[merged@meta.data$orig.ident=="F3",]$sampleType<-"MDD"
merged@meta.data[merged@meta.data$orig.ident=="F4",]$sampleType<-"MDD"
merged@meta.data[merged@meta.data$orig.ident=="F5",]$sampleType<-"MDD"
merged@meta.data[merged@meta.data$orig.ident=="F6",]$sampleType<-"MDD"
merged@meta.data[merged@meta.data$orig.ident=="F8",]$sampleType<-"MDD"
merged@meta.data[merged@meta.data$orig.ident=="F9",]$sampleType<-"MDD"
merged@meta.data[merged@meta.data$orig.ident=="M1",]$sampleType<-"MDD"
merged@meta.data[merged@meta.data$orig.ident=="M4",]$sampleType<-"MDD"
merged@meta.data[merged@meta.data$orig.ident=="M5",]$sampleType<-"MDD"
merged@meta.data[merged@meta.data$orig.ident=="M6",]$sampleType<-"MDD"
merged@meta.data[merged@meta.data$orig.ident=="M8",]$sampleType<-"MDD"
merged@meta.data[merged@meta.data$orig.ident=="M10",]$sampleType<-"MDD"
merged@meta.data[merged@meta.data$orig.ident=="M11",]$sampleType<-"MDD"
merged@meta.data[merged@meta.data$orig.ident=="M14",]$sampleType<-"MDD"
merged@meta.data[merged@meta.data$orig.ident=="M17",]$sampleType<-"MDD"
merged@meta.data[merged@meta.data$orig.ident=="M18",]$sampleType<-"MDD"
merged@meta.data[merged@meta.data$orig.ident=="M23",]$sampleType<-"MDD"
merged@meta.data[merged@meta.data$orig.ident=="M26",]$sampleType<-"MDD"
merged@meta.data[merged@meta.data$orig.ident=="M28",]$sampleType<-"MDD"
merged@meta.data[merged@meta.data$orig.ident=="M30",]$sampleType<-"MDD"
merged@meta.data[merged@meta.data$orig.ident=="M32",]$sampleType<-"MDD"
merged@meta.data[merged@meta.data$orig.ident=="M33",]$sampleType<-"MDD"
merged@meta.data[merged@meta.data$orig.ident=="M34",]$sampleType<-"MDD"
DefaultAssay(merged)<-"RNA"
merged<-JoinLayers(merged)
#annotaion for cell annotation
#all neuron (Ex_Neuron and In_Neuron)
new.cluster.ids <- c("celltype","Neurons","Neurons","Neurons","Neurons","Neurons","Neurons","Astrocyte","Astrocyte","Neurons","Neurons","Neurons","Neurons","Neurons","Neurons","Astrocyte","Astrocyte","Neurons","Neurons","Neurons","Neurons","Neurons","Neurons","Neurons","Neurons","Endothelial_cells","Endothelial_cells","Monocyte","Monocyte","Neurons","Neurons","Neurons","Neurons","Neurons","Neurons","Neurons","Neurons")
names(new.cluster.ids) <- levels(merged)
merged <- RenameIdents(merged, new.cluster.ids)
merged@meta.data$cell_anno <- Idents(merged)
#get cluster1.expr matrix
##
cluster1.list<-c("DENND1A","TOX","ERBB4","SNRK","PLCL2","NKAIN2")
##
cluster1.matrix<-merged@assays$RNA$data[cluster1.list,]
cluster1.matrix.cv<-t(as.data.frame(cluster1.matrix))
merged@meta.data<-cbind(merged@meta.data,cluster1.matrix.cv)
Idents(merged)<-merged@meta.data[,"cell_anno"]
neuron<-subset(merged,idents="Neurons")
#save metadata
saveRDS(neuron@meta.data,"SPD.cluster1.metadata.Neuron.rds")
#save only neuron seurat object
saveRDS(neuron,"Neuron.rds")


##
Rscript 2024-9-30.get_female_and_allPatient.r
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>figure plot analysis (visualization)<<<<<<<<<
cat >Dennd1a.visualization.r
meta<-readRDS("SPD.cluster1.metadata.Neuron.rds")
meta$sampleType2<-meta$orig.ident
meta[grepl("WT",meta$sampleType),]$sampleType2<-"WT"
meta$sampleType2<-factor(meta$sampleType2,levels = c('F1','F2','F3','F4','F5','F6','F8','F9',
                                                     'F11','F12','F14','F15','F16','F17','F18',
                                                     'F19','F20','F25','F27','F28','M1','M4',
                                                     'M5','M6','M8','M10','M11','M14','M17',
                                                     'M18','M23','M26','M28','M30','M32','M33','M34','WT'))

meta2<-meta[grepl("F",meta$orig.ident),]
table(meta2[grepl("WT",meta2$sampleType2),]$orig.ident)

table(meta$sampleType2)

#only F9,F14£¬F15 and Ctrl
meta$subSample<-"no"
meta[meta$orig.ident=="F9",]$subSample<-"yes"
meta[meta$orig.ident=="F14",]$subSample<-"yes"
meta[meta$orig.ident=="F15",]$subSample<-"yes"
meta[grepl("F",meta$orig.ident) & meta$sampleType=="WT",]$subSample<-"yes"
library(ggplot2)
#select Dennd1a as example
female_each_DENND1A<-ggplot(meta[grepl("F",meta$orig.ident),], aes(x=sampleType2, y=DENND1A)) + 
  geom_boxplot() + 
  stat_summary(fun=mean, geom="line", aes(group=1),col="red")  + 
  stat_summary(fun=mean, geom="point",col="red")+
  stat_summary(fun=median, geom="line", aes(group=1),col="green")  + 
  stat_summary(fun=median, geom="point",col="green")+
  geom_hline(yintercept=summary(meta[grepl("WT",meta$sampleType2) & grepl("F",meta$orig.ident),]$DENND1A)[4], color="#FCCDE5", linetype="dashed")+
  geom_hline(yintercept=summary(meta[grepl("WT",meta$sampleType2) & grepl("F",meta$orig.ident),]$DENND1A)[3], color="#CCEBC5", linetype="dashed")

female_each_DENND1A_sub<-ggplot(meta[meta$subSample=="yes",], aes(x=sampleType2, y=DENND1A)) + 
  geom_boxplot() + 
  stat_summary(fun=mean, geom="line", aes(group=1),col="red")  + 
  stat_summary(fun=mean, geom="point",col="red")+
  stat_summary(fun=median, geom="line", aes(group=1),col="green")  + 
  stat_summary(fun=median, geom="point",col="green")+
  geom_hline(yintercept=summary(meta[grepl("WT",meta$sampleType2) & grepl("F",meta$orig.ident),]$DENND1A)[4], color="#FCCDE5", linetype="dashed")+
  geom_hline(yintercept=summary(meta[grepl("WT",meta$sampleType2) & grepl("F",meta$orig.ident),]$DENND1A)[3], color="#CCEBC5", linetype="dashed")

pdf("female.neuron.each.cluster1_genes.check.pdf",width = 9.24,height = 3.99)
female_each_DENND1A
dev.off()
pdf("female.neuron.each.cluster1_genes.sub.check.pdf",width = 2.71,height = 2.74)
female_each_DENND1A_sub
dev.off()

#########
Rscript Dennd1a.visualization.r

#>>>>>>>>>>>>>>>>>>>>>>>GSEA analysis for Female sample
#get DEGs first, this part using F9 as example
(R4.3.1) xuzhzh@xdzpefect:/mnt/disk12_part2/zls_project/publicData_SPD/resolution0.1$ cat get_F1-vs-Female_wt.DEGs.r
library(Seurat)
merged<-readRDS("Neuron.rds")
head(merged@meta.data)
DefaultAssay(merged)
Idents(merged)<-merged@meta.data[,"orig.ident"]
table(merged@meta.data$orig.ident)
#female-vs-wt
de <- FindMarkers(object = merged, ident.1 = "F9", ident.2 = c("F7","F10","F13","F21","F22","F23","F24","F26","F29","F30","F31","F32","F33","F34","F35","F36","F37","F38"), logfc.threshold = 0.15, min.pct = 0.25, test.use = "LR")
head(de)
write.table(de,"F9-vs-Female_wt.DEGs.txt",quote=F)

#Running
Rscript  get_F9-vs-Female_wt.DEGs.r

#GSEA analysis
rm(list=ls())
setwd("d:/data_analysis/zls_project/0SPD_dual-sg_pertubseq_project/0publicData/")

#2024-10-29 GSEA analysis
rm(list = ls())
library(ReactomePA)
#BiocManager::install("ReactomePA") 
library(tidyverse)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
#BiocManager::install("clusterProfiler")
library(biomaRt)
#BiocManager::install("biomaRt") 
library(enrichplot)
library(clusterProfiler)
#updateEnsDb(species = "Homo sapiens")

#F9
data<-read.table("GSEA_analysis/F9-vs-Female_wt.DEGs.txt",header=T)
data$geneSymbol<-row.names(data)
ref2<-read.table("GSEA_analysis/all_geneSymbol.convert.txt",header=T,sep="\t")
colnames(ref2)<-c("geneSymbol","ENTREZID","species","gene_full_name")
data<-merge(data,ref2,by="geneSymbol")
de<-data
aaa<-data.frame(ENTREZID=de$ENTREZID,logFC=de$avg_log2FC)
aaa<-aaa[order(-aaa$logFC),]
geneList=aaa[,2]
names(geneList)=as.character(aaa[,1])
#GSEA-KEGG
KEGG_gseresult<-gseKEGG(geneList,nPerm=1000,minGSSize = 10,maxGSSize = 1000,pvalueCutoff = 1)
#Save outcome
write.table(KEGG_gseresult,file="GSEA_analysis/F9_allDEGs_KEGG_gsearesult.csv",sep=",",row.names = TRUE)
saveRDS(KEGG_gseresult,"GSEA_analysis/F9.allDEGs_KEGG_gseresult.rds")
KEGG_gseresult<-readRDS("GSEA_analysis/F9.allDEGs_KEGG_gseresult.rds")

pdf("GSEA_analysis/F9.allDEGs_KEGG_GSEAplot2.pdf")
gseaplot2(KEGG_gseresult,1,pvalue_table = TRUE)
ridgeplot(KEGG_gseresult,27)
ridgeplot(KEGG_gseresult,10)
dev.off()

sessionInfo()

