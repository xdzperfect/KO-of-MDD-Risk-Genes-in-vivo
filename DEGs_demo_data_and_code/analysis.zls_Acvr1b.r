library(Seurat)
ind<-readRDS("zls_Acvr1b.rds")
ind <- subset(ind,celltype %in% c("In_neuron","Ex_neuron"))
DefaultAssay(ind) <- "RNA"
Idents(ind)<-ind@meta.data[,"per_gene"]
de <- FindMarkers(object = ind,
                        ident.1 = "zls_Acvr1b",
                        ident.2 = "safe_H",
logfc.threshold = 0.15,
                        min.pct = 0.25,
                        test.use = "LR")
write.table(de,"zls_Acvr1b.de.findmarker.txt",quote=F)
library(Seurat)
de<-read.table("zls_Acvr1b.de.findmarker.txt",header=T)

genes<-rownames(de[de$p_val_adj<0.05,])
if(length(genes)<3){
print ("zls_Acvr1b didn't lead to a significant transcriptional phenotype")
}else{
library(dplyr)
if(length(genes)>50){
genes<-de %>% top_n(-50,p_val_adj) %>% rownames(.)
}
exp<-data.frame(t(as.data.frame(ind[["RNA"]]$data[genes,])))
exp$cell_name<-rownames(exp)
ind@meta.data$cell_name<-rownames(ind@meta.data)
meta<-ind@meta.data[,c("cell_name","per_gene")]
names(meta)<-c("cell_name","per_gene")
exp$cell_name<-rownames(exp)
exp<-merge(exp,meta,by="cell_name")
rownames(exp)<-exp$cell_name
exp<-exp[,-1]
library(MASS)
model<-MASS::lda(factor(per_gene)~.,data=exp)
pred<-predict(model,newdata=exp[,!names(exp) %in% "per_gene"])
post<-data.frame(pred$class)
exp<-cbind(exp,post)
exp$new_label<-"safe_H"
exp$new_label[exp$per_gene == "zls_Acvr1b" & exp$pred.class=="zls_Acvr1b"]<-paste("zls_Acvr1b","OE",sep="_")
exp$new_label[exp$per_gene == "zls_Acvr1b" & exp$pred.class=="safe_H"]<-paste("zls_Acvr1b","NP",sep="_")
exp$cell_name<-rownames(exp)
av<-nrow(exp[exp$new_label %in%  paste("zls_Acvr1b","OE",sep="_"),])/nrow(exp[exp$per_gene %in% "zls_Acvr1b",])*100
print(paste(as.character(round(av)),"% of cells show a significant per_geneation for ","zls_Acvr1b",sep=""))
ind$new_label<-"hold"
ind$new_label[ind@meta.data[,"per_gene"] %in% "safe_H"]<-"safe_H"
ind$new_label[ind@meta.data[,"cell_name"] %in% exp$cell_name[exp$new_label %in% paste("zls_Acvr1b","OE",sep="_")]]<-paste("zls_Acvr1b","OE",sep="_")
ind$new_label[ind@meta.data[,"cell_name"] %in% exp$cell_name[exp$new_label %in% paste("zls_Acvr1b","NP",sep="_")]]<-paste("zls_Acvr1b","NP",sep="_")
saveRDS(ind,"zls_Acvr1b.add.new_label.rds")
}
library(Seurat)
library(dplyr)
library("edgeR")
seuset<-readRDS("zls_Acvr1b.add.new_label.rds")
lda_after<-as.data.frame(table(seuset@meta.data$new_label))
 write.table(lda_after,"zls_Acvr1b.cellnumber.txt",quote=F)
if(dim(table(seuset@meta.data[seuset@meta.data$new_label=="safe_H",]$id))[1]>3){
avg_exp<-rowMeans(seuset[["RNA"]]$data[,colnames(seuset[["RNA"]])%in% rownames(seuset@meta.data[seuset$per_gene=="safe_H",])])
avg_exp<-avg_exp[avg_exp>0.25]
avg_exp<-data.frame(avg_exp)
avg_exp$gene<-rownames(avg_exp)
expr<-seuset[["RNA"]]$counts
meta<-seuset@meta.data
names(meta)[names(meta)=="orig.ident"]<-"sample"
names(meta)[names(meta)=="new_label"]<-"per_geneation"
mm<-model.matrix(~0+per_geneation:sample,data=meta)
mat_mm<-expr %*% mm
keep_genes<-rowSums(mat_mm>0)>0
mat_mm<-mat_mm[keep_genes,]%>% as.data.frame()
colnames(mat_mm)=gsub("per_geneation|sample","",colnames(mat_mm))
keep_samples=colSums(mat_mm)>0
mat_mm<-mat_mm[,keep_samples]
bulk<-mat_mm[avg_exp$gene,]
meta<-data.frame(pseudobulk=colnames(bulk),umi_count=colSums(bulk))
meta$per_gene<-gsub("\\:.*","",meta$pseudobulk)
meta<-within(meta,per_gene<-relevel(as.factor(per_gene),ref="safe_H"))
counts<-bulk
metadata<-meta
m<-metadata[metadata$per_gene %in% c("safe_H","zls_Acvr1b_OE"),]
counts<-counts[,names(counts) %in% m$pseudobulk]
design<-model.matrix(pseudobulk~per_gene,data=m)
design<-design[,colSums(design)!=0]
y<-DGEList(counts=counts,group=m$per_gene)%>% calcNormFactors(method='TMM')%>% estimateDisp(design)
fit=glmFit(y,design=design)
test=glmLRT(fit)
t1<-topTags(test,n=Inf)%>% as.data.frame()
t1$gene<-row.names(t1)
res=t1%>%mutate(de_family='pseudobulk',de_method="edgeR",de_type="LRT")
write.table(res,"DEGs.zls_Acvr1b.edgeR.txt",quote=F,row.names=F)
}else{
        print("pay attention: the control type of safe_H is less than 3!!!!!",end="\n");
        print("Cancel analysis!!!",end="\n");
}

