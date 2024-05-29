library(Seurat)
library(Matrix)


list <- read.table("datalist.txt",sep="\t",header=T,row.names=1)


for(i in 1:nrow(list)){
if(i == 1){
merge <- Read10X(list[i,7])
colnames(merge) <- paste(colnames(merge),"_",i,sep="")
}else{
data <- Read10X(list[i,7])
colnames(data) <- paste(colnames(data),"_",i,sep="")
merge <- cbind(merge,data)
}
}

merge <- CreateSeuratObject(counts=merge,project="combined",names.field = 2, names.delim = "_")
mito.gene <-  rownames(merge)[grep("MT-",rownames(merge))]
percent.mito <- Matrix::colSums(merge@assays$RNA[mito.gene, ])/Matrix::colSums(merge@assays$RNA)
AddMetaData(merge,metadata=percent.mito,col.name="percent.mito") -> merge

age <- rep(0,ncol(merge))
names(age) <- colnames(merge)
stage <- age
sex <- age

for(i in 1:nrow(list)){
      age[which(merge@meta.data$orig.ident == i)] <- list[i,2]
      stage[which(merge@meta.data$orig.ident == i)] <- list[i,3]
      sex[which(merge@meta.data$orig.ident == i)] <- list[i,5]
}

AddMetaData(merge,metadata=age,col.name="Age") -> merge
AddMetaData(merge,metadata=stage,col.name="Stage") -> merge
AddMetaData(merge,metadata=sex,col.name="Sex") -> merge


merge <- subset(merge,subset= nFeature_RNA > 100 & nFeature_RNA < 7000 & nCount_RNA > 500 & nCount_RNA < 20000 & percent.mito < 0.05)

merge.list <- SplitObject(object=merge,split.by="orig.ident")

for(i in 1:length(merge.list)){
merge.list[[i]] <- NormalizeData(merge.list[[i]],verbose=F)
merge.list[[i]] <- FindVariableFeatures(merge.list[[i]],selection.method="vst",nfeatures=2000,verbose=F)
}

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
save(merge.list,file="merge_list.dat")
save(merge.anchors,file="merge_anchors.dat")


