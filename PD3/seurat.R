library(Seurat)
library(Matrix)
list <- read.table("datalist.txt",sep="\t",header=T,row.names=1)
rownames(list) <- list[,1]
for(i in 1:nrow(list)){
 print(rownames(list)[i])
 file <- paste("map/",rownames(list)[i],"/outs/filtered_feature_bc_matrix.h5",sep="")
 data <- Read10X_h5(file)
 colnames(data) <- paste(colnames(data),"_",rownames(list)[i],sep="")
 if(i == 1){
 all <- data
 }else{
 all <- cbind(all,data)
 }
}
merge <- CreateSeuratObject(counts=all,project="combined",names.field = 2, names.delim = "_")
mito.gene <-  rownames(merge)[grep("MT-",rownames(merge))]
percent.mito <- Matrix::colSums(merge@assays$RNA[mito.gene, ])/Matrix::colSums(merge@assays$RNA)
AddMetaData(merge,metadata=percent.mito,col.name="percent.mito") -> merge

for(i in 1:ncol(list)){
      AddMetaData(merge,metadata=list[merge$orig.ident,i],col.name=colnames(list)[i]) -> merge
}

merge <- subset(merge,subset= nFeature_RNA > 100 & nFeature_RNA < 7000 & nCount_RNA > 500 & nCount_RNA < 30000 & percent.mito < 0.1)

merge.list <- SplitObject(object=merge,split.by="orig.ident")

for(i in 1:length(merge.list)){
merge.list[[i]] <- NormalizeData(merge.list[[i]],verbose=F)
merge.list[[i]] <- FindVariableFeatures(merge.list[[i]],selection.method="vst",nfeatures=2000,verbose=F)
}
save(merge.list,file="merge_list.dat")

merge.anchors <- FindIntegrationAnchors(merge.list,dims=1:20)
merge.integrated <- IntegrateData(anchorset = merge.anchors,dims=1:20)

library(ggplot2)
library(cowplot)

DefaultAssay(object = merge.integrated) <- "integrated"
merge.integrated <- ScaleData(object = merge.integrated, verbose = FALSE)
merge.integrated <- RunPCA(object = merge.integrated, npcs = 20, verbose = FALSE)
merge.integrated <- RunUMAP(object = merge.integrated, reduction = "pca", dims = 1:20)

merge.integrated <- FindNeighbors(merge.integrated,dims=1:20,reduction="pca")
merge.integrated <- FindClusters(merge.integrated)

save(merge.integrated,file="merge.dat")

markers <- FindAllMarkers(merge.integrated)
save(markers,file="markers.dat")

new_group <- rep(0,ncol(merge.integrated))
names(new_group) <- colnames(merge.integrated)
new_group[which(merge.integrated@active.ident == 18)] <- 1
new_group[which(merge.integrated@active.ident == 3)] <- 1
new_group[which(merge.integrated@active.ident == 4)] <- 1
new_group[which(merge.integrated@active.ident == 13)] <- 1
new_group[which(merge.integrated@active.ident == 17)] <- 1
new_group[which(merge.integrated@active.ident == 20)] <- 1
new_group[which(merge.integrated@active.ident == 19)] <- 1
new_group[which(merge.integrated@active.ident == 15)] <- 1
new_group[which(merge.integrated@active.ident == 10)] <- 2
new_group[which(merge.integrated@active.ident == 2)] <- 2
new_group[which(merge.integrated@active.ident == 0)] <- 2
new_group[which(merge.integrated@active.ident == 1)] <- 2
new_group[which(merge.integrated@active.ident == 9)] <- 3
new_group[which(merge.integrated@active.ident == 6)] <- 4
new_group[which(merge.integrated@active.ident == 16)] <- 4
new_group[which(merge.integrated@active.ident == 7)] <- 4
new_group[which(merge.integrated@active.ident == 8)] <- 4
new_group[which(merge.integrated@active.ident == 11)] <- 4
new_group[which(merge.integrated@active.ident == 21)] <- 4
new_group[which(merge.integrated@active.ident == 12)] <- 5
new_group[which(merge.integrated@active.ident == 14)] <- 5
new_group[which(merge.integrated@active.ident == 22)] <- 5
new_group[which(merge.integrated@active.ident == 5)] <- 6
new_group[which(merge.integrated@active.ident == 23)] <- 7
merge.integrated@meta.data$CellType <- factor(new_group)

