#get Query mapped cells
Rscript get_mapped_cell.r

#part0 get each perturb genes' seurat object
cat >get_rds.r

library(Seurat)
merged<-readRDS("sPD.1moi.query2ref_MDD.matched.Ex_neuron.rds")
perturb<-as.data.frame(table(merged@meta.data$id))
perturblist<-subset(perturb,!grepl("WBL|ctrl|SPD",perturb$Var1) & Freq>2)$Var1
merged@meta.data$per_gene<-merged@meta.data$id
merged@meta.data[grepl("WBL|SPD|ctrl",merged@meta.data$id),]$per_gene<-"control"
Idents(merged)<-merged@meta.data[,"per_gene"]
for(i in perturblist){
ind<-subset(merged,idents=c(i,"control"))
saveRDS(ind,paste0("Ex_neuron/",i,".rds"))
}

#running script
Rscript get_rds.r

#part1 get DEGs by using findmarker functions

cd /home/xuzz01/mv/zls_project/2024-4-1-dual-sg-SPD/single2025-6-12-seuratV5_query/script_neurons/part1_de.findmarker
(base) -bash-4.2$ cat pipeline.sh
#!/bin/bash
#SBATCH -J pipeline
#SBATCH -n 1
#SBATCH --cpus-per-task 5
#SBATCH --mem=10
#SBATCH -o pipeline_output_%J
#SBATCH -e pipeline_errput_%J
#SBATCH -p basic
cd /home/xuzz01/mv/zls_project/2024-4-1-dual-sg-SPD/single2025-6-12-seuratV5_query/script_neurons/part1_de.findmarker
Rscript pipeline.r
(base) -bash-4.2$ cat pipeline.r
library(Seurat)
setwd("/home/xuzz01/mv/zls_project/2024-4-1-dual-sg-SPD/single2025-6-12-seuratV5_query/Ex_neuron")
ind<-readRDS("pipeline.rds")
Idents(ind)<-ind@meta.data[,"treatment"]
de <- FindMarkers(object = ind,
                        ident.1 = "treatment",
                        ident.2 = "control",
logfc.threshold = 0.15,
                        min.pct = 0.25,
                        test.use = "LR")
write.table(de,"pipeline.de.findmarker.txt",quote=F)


#
ls ../../Ex_neuron/*rds|awk -F 'Ex_neuron/' '{print $2}'|sed 's/.rds//g'|while read i
do
cp pipeline.r "$i".r
cp pipeline.sh "$i".sh
sed -i "s/pipeline/$i/g" "$i".*
sbatch "$i".sh
done 


#part2 pseudo bulk analysis
cd /home/xuzz01/mv/zls_project/2024-4-1-dual-sg-SPD/single2025-6-12-seuratV5_query/script_neurons/pattern0_paper_method_findmarker2LDA2pseudoBulk

(base) -bash-4.2$ cat pipeline.r
library(Seurat)
setwd("/home/xuzz01/mv/zls_project/2024-4-1-dual-sg-SPD/single2025-6-12-seuratV5_query/Ex_neuron")
de<-read.table("pipeline.de.findmarker.txt",header=T)
ind<-readRDS("pipeline.rds")
ind@meta.data$per_gene<-gsub("control","safe_H",ind@meta.data$per_gene)
genes<-rownames(de[de$p_val_adj<0.05,])
if(length(genes)<5){print ("pipeline didn't lead to a significant transcriptional phenotype")}
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
exp$new_label[exp$per_gene == "pipeline" & exp$pred.class=="pipeline"]<-paste("pipeline","OE",sep="_")
exp$new_label[exp$per_gene == "pipeline" & exp$pred.class=="safe_H"]<-paste("pipeline","NP",sep="_")
exp$cell_name<-rownames(exp)
av<-nrow(exp[exp$new_label %in%  paste("pipeline","OE",sep="_"),])/nrow(exp[exp$per_gene %in% "pipeline",])*100
print(paste(as.character(round(av)),"% of cells show a significant per_geneation for ","pipeline",sep=""))
ind$new_label<-"hold"
ind$new_label[ind@meta.data[,"per_gene"] %in% "safe_H"]<-"safe_H"
ind$new_label[ind@meta.data[,"cell_name"] %in% exp$cell_name[exp$new_label %in% paste("pipeline","OE",sep="_")]]<-paste("pipeline","OE",sep="_")
ind$new_label[ind@meta.data[,"cell_name"] %in% exp$cell_name[exp$new_label %in% paste("pipeline","NP",sep="_")]]<-paste("pipeline","NP",sep="_")
saveRDS(ind,"pipeline.add.new_label.rds")
#part2
library(Seurat)
library(dplyr)
library("edgeR")
seuset<-readRDS("pipeline.add.new_label.rds")
lda_after<-as.data.frame(table(seuset@meta.data$new_label))
write.table(lda_after,"pipeline.cellnumber.txt",quote=F)
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
m<-metadata[metadata$per_gene %in% c("safe_H","pipeline_OE"),]
counts<-counts[,names(counts) %in% m$pseudobulk]
design<-model.matrix(pseudobulk~per_gene,data=m)
design<-design[,colSums(design)!=0]
y<-DGEList(counts=counts,group=m$per_gene)%>% calcNormFactors(method='TMM')%>% estimateDisp(design)
fit=glmFit(y,design=design)
test=glmLRT(fit)
t1<-topTags(test,n=Inf)%>% as.data.frame()
t1$gene<-row.names(t1)
res=t1%>%mutate(de_family='pseudobulk',de_method="edgeR",de_type="LRT")
write.table(res,"DEGs.pipeline.edgeR.txt",quote=F,row.names=F)

(base) -bash-4.2$ cat pipeline.sh
#!/bin/bash
#SBATCH -J pipeline
#SBATCH -n 1
#SBATCH --cpus-per-task 2
#SBATCH --mem=10
#SBATCH -o pipeline_output_%J
#SBATCH -e pipeline_errput_%J
#SBATCH -p basic
cd /home/xuzz01/mv/zls_project/2024-4-1-dual-sg-SPD/single2025-6-12-seuratV5_query/script_neurons/pattern0_paper_method_findmarker2LDA2pseudoBulk
Rscript pipeline.r

#running pipeline
ls -l ../../Ex_neuron/*.de.findmarker.txt|awk -F 'Ex_neuron/' '{print $2}'|sed 's/.de.findmarker.txt//g'|while read i
do
cp pipeline.r "$i".r
cp pipeline.sh "$i".sh
sed -i "s/pipeline/$i/g" "$i".*
sbatch "$i".sh
done 

#part3 pathway enrichment analysis using WebGestaltR and gprofiler2 R packages


