rm(list=ls())
setwd("d:/data_analysis/zls_project/0SPD_dual-sg_pertubseq_project/2025-Nature_genetics_respone_prepare/2025-6-3-specific_pathway/")


library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(magick)
library(dendsort)
library(dplyr)
library(tidyverse)
###
data<-read.csv("final.csv",header=T)
head(data[,1:10])
dim(data)
#GABAergic synapase
ref<-read.table("GABAergic_synapase.txt",header=T)
colnames(ref)<-"geneID"
colnames(ref)

head(ref)
data2<-merge(data,ref,by="geneID")
data2<-data2%>%column_to_rownames("geneID")
write.csv(data2,"GABAergic_synapase.data.csv",quote=F)
tmp<-as.data.frame(t(data2))
class(tmp)
tmp$sum<-apply(tmp,1,sum)
dim(tmp)
#filter sum==0
tmp<-tmp[tmp$sum!=0,]
dim(tmp)
data2<-as.data.frame(t(tmp[,1:89]))
class(data2)
cor<-cor(data2)
write.csv(cor,"GABAergic_synapase.cordata.csv")

#get order loci
cluster=hclust(dist(cor))

p0<-
  Heatmap(cor, 
          name = "Persion\ncorrelation", 
          col = colorRamp2(c(-1, 0, 1), c("#081d58", "#41b6c4", "yellow")), 
          show_row_names = T, 
          cluster_columns=cluster,
          show_column_names = T, 
          show_row_dend = T, 
          show_column_dend = T,
          cluster_rows=cluster,
          column_title = "SPD perturbation in GABAergic synapase")

p1<-
  Heatmap(cor, 
          name = "Persion\ncorrelation", 
          col = colorRamp2(c(0.4, 0.6, 0.8), c("#081d58", "#41b6c4", "yellow")), 
          show_row_names = T, 
          cluster_columns=cluster,
          show_column_names = T, 
          show_row_dend = T, 
          show_column_dend = T,
          cluster_rows=cluster,
          column_title = "SPD perturbation in GABAergic synapase")
p2<-
  Heatmap(cor, 
          name = "Persion\ncorrelation", 
          col = colorRamp2(c(0.4, 0.5, 0.6), c("#081d58", "#41b6c4", "yellow")), 
          show_row_names = T, 
          cluster_columns=cluster,
          show_column_names = T, 
          show_row_dend = T, 
          show_column_dend = T,
          cluster_rows=cluster,
          column_title = "SPD perturbation in GABAergic synapase")
p3<-
  Heatmap(cor, 
          name = "Persion\ncorrelation", 
          col = colorRamp2(c(0, 0.25, 0.5), c("#081d58", "#41b6c4", "yellow")), 
          show_row_names = T, 
          cluster_columns=cluster,
          show_column_names = T, 
          show_row_dend = T, 
          show_column_dend = T,
          cluster_rows=cluster,
          column_title = "SPD perturbation in GABAergic synapase")
pdf("GABAergic_synapase.pdf")
p0
p1
p2
p3
dev.off()
png("GABAergic_synapase0.png")
p0
dev.off()
png("GABAergic_synapase1.png")
p1
dev.off()
png("GABAergic_synapase2.png")
p2
dev.off()
png("GABAergic_synapase3.png")
p3
dev.off()

saveRDS(cor,"GABAergic_synapase.cor.rds")
saveRDS(cluster,"GABAergic_synapase.cor_cluster.rds")

#
#Dopaminergic synapse
ref<-read.table("Dopaminergic_synapse.txt",header=T)
colnames(ref)<-"geneID"
colnames(ref)
dim(ref)
ref<-unique(ref)
head(ref)
data2<-merge(data,ref,by="geneID")
data2<-data2%>%column_to_rownames("geneID")
write.csv(data2,"Dopaminergic_synapse.data.csv",quote=F)
tmp<-as.data.frame(t(data2))
class(tmp)
tmp$sum<-apply(tmp,1,sum)
dim(tmp)
#filter sum==0
tmp<-tmp[tmp$sum!=0,]
dim(tmp)
data2<-as.data.frame(t(tmp[,1:132]))
class(data2)
cor<-cor(data2)
write.csv(cor,"Dopaminergic_synapse.cordata.csv")

#get order loci
cluster=hclust(dist(cor))

p0<-
  Heatmap(cor, 
          name = "Persion\ncorrelation", 
          col = colorRamp2(c(-1, 0, 1), c("#081d58", "#41b6c4", "yellow")), 
          show_row_names = T, 
          cluster_columns=cluster,
          show_column_names = T, 
          show_row_dend = T, 
          show_column_dend = T,
          cluster_rows=cluster,
          column_title = "SPD perturbation in Dopaminergic synapse")

p1<-
  Heatmap(cor, 
          name = "Persion\ncorrelation", 
          col = colorRamp2(c(0.4, 0.6, 0.8), c("#081d58", "#41b6c4", "yellow")), 
          show_row_names = T, 
          cluster_columns=cluster,
          show_column_names = T, 
          show_row_dend = T, 
          show_column_dend = T,
          cluster_rows=cluster,
          column_title = "SPD perturbation in Dopaminergic synapse")
p2<-
  Heatmap(cor, 
          name = "Persion\ncorrelation", 
          col = colorRamp2(c(0.4, 0.5, 0.6), c("#081d58", "#41b6c4", "yellow")), 
          show_row_names = T, 
          cluster_columns=cluster,
          show_column_names = T, 
          show_row_dend = T, 
          show_column_dend = T,
          cluster_rows=cluster,
          column_title = "SPD perturbation in Dopaminergic synapse")
p3<-
  Heatmap(cor, 
          name = "Persion\ncorrelation", 
          col = colorRamp2(c(0, 0.25, 0.5), c("#081d58", "#41b6c4", "yellow")), 
          show_row_names = T, 
          cluster_columns=cluster,
          show_column_names = T, 
          show_row_dend = T, 
          show_column_dend = T,
          cluster_rows=cluster,
          column_title = "SPD perturbation in Dopaminergic synapse")
pdf("Dopaminergic_synapse.pdf")
p0
p1
p2
p3
dev.off()
png("Dopaminergic_synapse0.png")
p0
dev.off()
png("Dopaminergic_synapse1.png")
p1
dev.off()
png("Dopaminergic_synapse2.png")
p2
dev.off()
png("Dopaminergic_synapse3.png")
p3
dev.off()

saveRDS(cor,"Dopaminergic_synapse.cor.rds")
saveRDS(cluster,"Dopaminergic_synapse.cor_cluster.rds")

#Glutamatergic synapsee
ref<-read.table("Glutamatergic_synapse.txt",header=T)
colnames(ref)<-"geneID"
colnames(ref)

head(ref)
data2<-merge(data,ref,by="geneID")
data2<-data2%>%column_to_rownames("geneID")
write.csv(data2,"Glutamatergic_synapse.data.csv",quote=F)
tmp<-as.data.frame(t(data2))
class(tmp)
tmp$sum<-apply(tmp,1,sum)
dim(tmp)
#filter sum==0
tmp<-tmp[tmp$sum!=0,]
dim(tmp)
data2<-as.data.frame(t(tmp[,1:115]))
class(data2)
cor<-cor(data2)
write.csv(cor,"Glutamatergic_synapse.cordata.csv")

#get order loci
cluster=hclust(dist(cor))

p0<-
  Heatmap(cor, 
          name = "Persion\ncorrelation", 
          col = colorRamp2(c(-1, 0, 1), c("#081d58", "#41b6c4", "yellow")), 
          show_row_names = T, 
          cluster_columns=cluster,
          show_column_names = T, 
          show_row_dend = T, 
          show_column_dend = T,
          cluster_rows=cluster,
          column_title = "SPD perturbation in Glutamatergic synapse")

p1<-
  Heatmap(cor, 
          name = "Persion\ncorrelation", 
          col = colorRamp2(c(0.4, 0.6, 0.8), c("#081d58", "#41b6c4", "yellow")), 
          show_row_names = T, 
          cluster_columns=cluster,
          show_column_names = T, 
          show_row_dend = T, 
          show_column_dend = T,
          cluster_rows=cluster,
          column_title = "SPD perturbation in Glutamatergic synapse")
p2<-
  Heatmap(cor, 
          name = "Persion\ncorrelation", 
          col = colorRamp2(c(0.4, 0.5, 0.6), c("#081d58", "#41b6c4", "yellow")), 
          show_row_names = T, 
          cluster_columns=cluster,
          show_column_names = T, 
          show_row_dend = T, 
          show_column_dend = T,
          cluster_rows=cluster,
          column_title = "SPD perturbation in Glutamatergic synapse")
p3<-
  Heatmap(cor, 
          name = "Persion\ncorrelation", 
          col = colorRamp2(c(0, 0.25, 0.5), c("#081d58", "#41b6c4", "yellow")), 
          show_row_names = T, 
          cluster_columns=cluster,
          show_column_names = T, 
          show_row_dend = T, 
          show_column_dend = T,
          cluster_rows=cluster,
          column_title = "SPD perturbation in Glutamatergic synapse")
pdf("Glutamatergic_synapse.pdf")
p0
p1
p2
p3
dev.off()
png("Glutamatergic_synapse0.png")
p0
dev.off()
png("Glutamatergic_synapse1.png")
p1
dev.off()
png("Glutamatergic_synapse2.png")
p2
dev.off()
png("Glutamatergic_synapse3.png")
p3
dev.off()

saveRDS(cor,"Glutamatergic_synapse.cor.rds")
saveRDS(cluster,"Glutamatergic_synapse.cor_cluster.rds")

#Cholinergic synapse
ref<-read.table("Cholinergic_synapse.txt",header=T)
colnames(ref)<-"geneID"
colnames(ref)

head(ref)
data2<-merge(data,ref,by="geneID")
data2<-data2%>%column_to_rownames("geneID")
write.csv(data2,"Cholinergic_synapse.data.csv",quote=F)
tmp<-as.data.frame(t(data2))
class(tmp)
tmp$sum<-apply(tmp,1,sum)
dim(tmp)
#filter sum==0
tmp<-tmp[tmp$sum!=0,]
dim(tmp)
data2<-as.data.frame(t(tmp[,1:112]))
class(data2)
cor<-cor(data2)
write.csv(cor,"Cholinergic_synapse.cordata.csv")

#get order loci
cluster=hclust(dist(cor))

p0<-
  Heatmap(cor, 
          name = "Persion\ncorrelation", 
          col = colorRamp2(c(-1, 0, 1), c("#081d58", "#41b6c4", "yellow")), 
          show_row_names = T, 
          cluster_columns=cluster,
          show_column_names = T, 
          show_row_dend = T, 
          show_column_dend = T,
          cluster_rows=cluster,
          column_title = "SPD perturbation in Cholinergic synapse")

p1<-
  Heatmap(cor, 
          name = "Persion\ncorrelation", 
          col = colorRamp2(c(0.4, 0.6, 0.8), c("#081d58", "#41b6c4", "yellow")), 
          show_row_names = T, 
          cluster_columns=cluster,
          show_column_names = T, 
          show_row_dend = T, 
          show_column_dend = T,
          cluster_rows=cluster,
          column_title = "SPD perturbation in Cholinergic synapse")
p2<-
  Heatmap(cor, 
          name = "Persion\ncorrelation", 
          col = colorRamp2(c(0.4, 0.5, 0.6), c("#081d58", "#41b6c4", "yellow")), 
          show_row_names = T, 
          cluster_columns=cluster,
          show_column_names = T, 
          show_row_dend = T, 
          show_column_dend = T,
          cluster_rows=cluster,
          column_title = "SPD perturbation in Cholinergic synapse")
p3<-
  Heatmap(cor, 
          name = "Persion\ncorrelation", 
          col = colorRamp2(c(0, 0.25, 0.5), c("#081d58", "#41b6c4", "yellow")), 
          show_row_names = T, 
          cluster_columns=cluster,
          show_column_names = T, 
          show_row_dend = T, 
          show_column_dend = T,
          cluster_rows=cluster,
          column_title = "SPD perturbation in Cholinergic synapse")
pdf("Cholinergic_synapse.pdf")
p0
p1
p2
p3
dev.off()
png("Cholinergic_synapse0.png")
p0
dev.off()
png("Cholinergic_synapse1.png")
p1
dev.off()
png("Cholinergic_synapse2.png")
p2
dev.off()
png("Cholinergic_synapse3.png")
p3
dev.off()

saveRDS(cor,"Cholinergic_synapse.cor.rds")
saveRDS(cluster,"Cholinergic_synapse.cor_cluster.rds")

#Synaptic vesicle cycle
ref<-read.table("Synaptic_vesicle_cycle.txt",header=T)
colnames(ref)<-"geneID"
colnames(ref)

head(ref)
data2<-merge(data,ref,by="geneID")
data2<-data2%>%column_to_rownames("geneID")
write.csv(data2,"Synaptic_vesicle_cycle.data.csv",quote=F)
tmp<-as.data.frame(t(data2))
class(tmp)
tmp$sum<-apply(tmp,1,sum)
dim(tmp)
#filter sum==0
tmp<-tmp[tmp$sum!=0,]
dim(tmp)
data2<-as.data.frame(t(tmp[,1:77]))
class(data2)
cor<-cor(data2)
write.csv(cor,"Synaptic_vesicle_cycle.cordata.csv")

#get order loci
cluster=hclust(dist(cor))

p0<-
  Heatmap(cor, 
          name = "Persion\ncorrelation", 
          col = colorRamp2(c(-1, 0, 1), c("#081d58", "#41b6c4", "yellow")), 
          show_row_names = T, 
          cluster_columns=cluster,
          show_column_names = T, 
          show_row_dend = T, 
          show_column_dend = T,
          cluster_rows=cluster,
          column_title = "SPD perturbation in Synaptic vesicle cycle")

p1<-
  Heatmap(cor, 
          name = "Persion\ncorrelation", 
          col = colorRamp2(c(0.4, 0.6, 0.8), c("#081d58", "#41b6c4", "yellow")), 
          show_row_names = T, 
          cluster_columns=cluster,
          show_column_names = T, 
          show_row_dend = T, 
          show_column_dend = T,
          cluster_rows=cluster,
          column_title = "SPD perturbation in Synaptic vesicle cycle")
p2<-
  Heatmap(cor, 
          name = "Persion\ncorrelation", 
          col = colorRamp2(c(0.4, 0.5, 0.6), c("#081d58", "#41b6c4", "yellow")), 
          show_row_names = T, 
          cluster_columns=cluster,
          show_column_names = T, 
          show_row_dend = T, 
          show_column_dend = T,
          cluster_rows=cluster,
          column_title = "SPD perturbation in Synaptic vesicle cycle")
p3<-
  Heatmap(cor, 
          name = "Persion\ncorrelation", 
          col = colorRamp2(c(0, 0.25, 0.5), c("#081d58", "#41b6c4", "yellow")), 
          show_row_names = T, 
          cluster_columns=cluster,
          show_column_names = T, 
          show_row_dend = T, 
          show_column_dend = T,
          cluster_rows=cluster,
          column_title = "SPD perturbation in Synaptic vesicle cycle")
pdf("Synaptic_vesicle_cycle.pdf")
p0
p1
p2
p3
dev.off()
png("Synaptic_vesicle_cycle0.png")
p0
dev.off()
png("Synaptic_vesicle_cycle1.png")
p1
dev.off()
png("Synaptic_vesicle_cycle2.png")
p2
dev.off()
png("Synaptic_vesicle_cycle3.png")
p3
dev.off()

saveRDS(cor,"Synaptic_vesicle_cycle.cor.rds")
saveRDS(cluster,"Synaptic_vesicle_cycle.cor_cluster.rds")

#Endocytosis
ref<-read.table("Endocytosis.txt",header=T)
colnames(ref)<-"geneID"
colnames(ref)
ref<-unique(ref)
head(ref)
data2<-merge(data,ref,by="geneID")
data2$geneID<-gsub("-","_",data2$geneID)

data2<-data2%>%column_to_rownames("geneID")
write.csv(data2,"Endocytosis.data.csv",quote=F)
tmp<-as.data.frame(t(data2))
class(tmp)
tmp$sum<-apply(tmp,1,sum)
dim(tmp)
#filter sum==0
tmp<-tmp[tmp$sum!=0,]
dim(tmp)
data2<-as.data.frame(t(tmp[,1:259]))
class(data2)
cor<-cor(data2)
write.csv(cor,"Endocytosis.cordata.csv")

#get order loci
cluster=hclust(dist(cor))

p0<-
  Heatmap(cor, 
          name = "Persion\ncorrelation", 
          col = colorRamp2(c(-1, 0, 1), c("#081d58", "#41b6c4", "yellow")), 
          show_row_names = T, 
          cluster_columns=cluster,
          show_column_names = T, 
          show_row_dend = T, 
          show_column_dend = T,
          cluster_rows=cluster,
          column_title = "SPD perturbation in Endocytosis")

p1<-
  Heatmap(cor, 
          name = "Persion\ncorrelation", 
          col = colorRamp2(c(0.4, 0.6, 0.8), c("#081d58", "#41b6c4", "yellow")), 
          show_row_names = T, 
          cluster_columns=cluster,
          show_column_names = T, 
          show_row_dend = T, 
          show_column_dend = T,
          cluster_rows=cluster,
          column_title = "SPD perturbation in Endocytosis")
p2<-
  Heatmap(cor, 
          name = "Persion\ncorrelation", 
          col = colorRamp2(c(0.4, 0.5, 0.6), c("#081d58", "#41b6c4", "yellow")), 
          show_row_names = T, 
          cluster_columns=cluster,
          show_column_names = T, 
          show_row_dend = T, 
          show_column_dend = T,
          cluster_rows=cluster,
          column_title = "SPD perturbation in Endocytosis")
p3<-
  Heatmap(cor, 
          name = "Persion\ncorrelation", 
          col = colorRamp2(c(0, 0.25, 0.5), c("#081d58", "#41b6c4", "yellow")), 
          show_row_names = T, 
          cluster_columns=cluster,
          show_column_names = T, 
          show_row_dend = T, 
          show_column_dend = T,
          cluster_rows=cluster,
          column_title = "SPD perturbation in Endocytosis")
pdf("Endocytosis.pdf")
p0
p1
p2
p3
dev.off()
png("Endocytosis0.png")
p0
dev.off()
png("Endocytosis1.png")
p1
dev.off()
png("Endocytosis2.png")
p2
dev.off()
png("Endocytosis3.png")
p3
dev.off()

saveRDS(cor,"Endocytosis.cor.rds")
saveRDS(cluster,"Endocytosis.cor_cluster.rds")

#Axon guidance
ref<-read.table("Axon_guidance.txt",header=T)
colnames(ref)<-"geneID"
colnames(ref)
ref<-unique(ref)
head(ref)
data2<-merge(data,ref,by="geneID")
data2$geneID<-gsub("-","_",data2$geneID)

data2<-data2%>%column_to_rownames("geneID")
write.csv(data2,"Axon_guidance.data.csv",quote=F)
tmp<-as.data.frame(t(data2))
class(tmp)
tmp$sum<-apply(tmp,1,sum)
dim(tmp)
#filter sum==0
tmp<-tmp[tmp$sum!=0,]
dim(tmp)
data2<-as.data.frame(t(tmp[,1:181]))
class(data2)
cor<-cor(data2)
write.csv(cor,"Axon_guidance.cordata.csv")

#get order loci
cluster=hclust(dist(cor))

p0<-
  Heatmap(cor, 
          name = "Persion\ncorrelation", 
          col = colorRamp2(c(-1, 0, 1), c("#081d58", "#41b6c4", "yellow")), 
          show_row_names = T, 
          cluster_columns=cluster,
          show_column_names = T, 
          show_row_dend = T, 
          show_column_dend = T,
          cluster_rows=cluster,
          column_title = "SPD perturbation in Axon guidance")

p1<-
  Heatmap(cor, 
          name = "Persion\ncorrelation", 
          col = colorRamp2(c(0.4, 0.6, 0.8), c("#081d58", "#41b6c4", "yellow")), 
          show_row_names = T, 
          cluster_columns=cluster,
          show_column_names = T, 
          show_row_dend = T, 
          show_column_dend = T,
          cluster_rows=cluster,
          column_title = "SPD perturbation in Axon guidance")
p2<-
  Heatmap(cor, 
          name = "Persion\ncorrelation", 
          col = colorRamp2(c(0.4, 0.5, 0.6), c("#081d58", "#41b6c4", "yellow")), 
          show_row_names = T, 
          cluster_columns=cluster,
          show_column_names = T, 
          show_row_dend = T, 
          show_column_dend = T,
          cluster_rows=cluster,
          column_title = "SPD perturbation in Axon guidance")
p3<-
  Heatmap(cor, 
          name = "Persion\ncorrelation", 
          col = colorRamp2(c(0, 0.25, 0.5), c("#081d58", "#41b6c4", "yellow")), 
          show_row_names = T, 
          cluster_columns=cluster,
          show_column_names = T, 
          show_row_dend = T, 
          show_column_dend = T,
          cluster_rows=cluster,
          column_title = "SPD perturbation in Axon guidance")
pdf("Axon_guidance.pdf")
p0
p1
p2
p3
dev.off()
png("Axon_guidance0.png")
p0
dev.off()
png("Axon_guidance1.png")
p1
dev.off()
png("Axon_guidance2.png")
p2
dev.off()
png("Axon_guidance3.png")
p3
dev.off()

saveRDS(cor,"Axon_guidance.cor.rds")
saveRDS(cluster,"Axon_guidance.cor_cluster.rds")

#Calcium signaling pathway
ref<-read.table("Calcium_signaling_pathway.txt",header=T)
colnames(ref)<-"geneID"
colnames(ref)
ref<-unique(ref)
head(ref)
data2<-merge(data,ref,by="geneID")
data2$geneID<-gsub("-","_",data2$geneID)

data2<-data2%>%column_to_rownames("geneID")
write.csv(data2,"Calcium_signaling_pathway.data.csv",quote=F)
tmp<-as.data.frame(t(data2))
class(tmp)
tmp$sum<-apply(tmp,1,sum)
dim(tmp)
#filter sum==0
tmp<-tmp[tmp$sum!=0,]
dim(tmp)
data2<-as.data.frame(t(tmp[,1:253]))
class(data2)
cor<-cor(data2)
write.csv(cor,"Calcium_signaling_pathway.cordata.csv")

#get order loci
cluster=hclust(dist(cor))

p0<-
  Heatmap(cor, 
          name = "Persion\ncorrelation", 
          col = colorRamp2(c(-1, 0, 1), c("#081d58", "#41b6c4", "yellow")), 
          show_row_names = T, 
          cluster_columns=cluster,
          show_column_names = T, 
          show_row_dend = T, 
          show_column_dend = T,
          cluster_rows=cluster,
          column_title = "SPD perturbation in Calcium signaling pathway")

p1<-
  Heatmap(cor, 
          name = "Persion\ncorrelation", 
          col = colorRamp2(c(0.4, 0.6, 0.8), c("#081d58", "#41b6c4", "yellow")), 
          show_row_names = T, 
          cluster_columns=cluster,
          show_column_names = T, 
          show_row_dend = T, 
          show_column_dend = T,
          cluster_rows=cluster,
          column_title = "SPD perturbation in Calcium signaling pathway")
p2<-
  Heatmap(cor, 
          name = "Persion\ncorrelation", 
          col = colorRamp2(c(0.4, 0.5, 0.6), c("#081d58", "#41b6c4", "yellow")), 
          show_row_names = T, 
          cluster_columns=cluster,
          show_column_names = T, 
          show_row_dend = T, 
          show_column_dend = T,
          cluster_rows=cluster,
          column_title = "SPD perturbation in Calcium signaling pathway")
p3<-
  Heatmap(cor, 
          name = "Persion\ncorrelation", 
          col = colorRamp2(c(0, 0.25, 0.5), c("#081d58", "#41b6c4", "yellow")), 
          show_row_names = T, 
          cluster_columns=cluster,
          show_column_names = T, 
          show_row_dend = T, 
          show_column_dend = T,
          cluster_rows=cluster,
          column_title = "SPD perturbation in Calcium signaling pathway")
pdf("Calcium_signaling_pathway.pdf")
p0
p1
p2
p3
dev.off()
png("Calcium_signaling_pathway0.png")
p0
dev.off()
png("Calcium_signaling_pathway1.png")
p1
dev.off()
png("Calcium_signaling_pathway2.png")
p2
dev.off()
png("Calcium_signaling_pathway3.png")
p3
dev.off()

saveRDS(cor,"Axon_guidance.cor.rds")
saveRDS(cluster,"Axon_guidance.cor_cluster.rds")
