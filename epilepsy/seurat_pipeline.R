library(Seurat)
library(Matrix)

merge <- Read10X_h5("merge/outs/count/filtered_feature_bc_matrix.h5")
merge <- CreateSeuratObject(counts=merge,project="combined",names.field = 2, names.delim = "-")
mito.gene <-  rownames(merge)[grep("mt-",rownames(merge))]
percent.mito <- Matrix::colSums(merge@assays$RNA[mito.gene, ])/Matrix::colSums(merge@assays$RNA)
AddMetaData(merge,metadata=percent.mito,col.name="percent.mito") -> merge

type <- rep("control",ncol(merge))
type[which(merge@meta.data$orig.ident == 11)] <- "epilepsy"
type[which(merge@meta.data$orig.ident == 14)] <- "epilepsy"
type[which(merge@meta.data$orig.ident == 16)] <- "epilepsy"
type[which(merge@meta.data$orig.ident == 17)] <- "epilepsy"
type[which(merge@meta.data$orig.ident == 20)] <- "epilepsy"
type[which(merge@meta.data$orig.ident == 21)] <- "epilepsy"
type[which(merge@meta.data$orig.ident == 22)] <- "epilepsy"
type[which(merge@meta.data$orig.ident == 23)] <- "epilepsy"
AddMetaData(merge,metadata=type,col.name="Type") -> merge

age <- rep(19,ncol(merge))
age[which(merge@meta.data$orig.ident == 1)] <- 15
age[which(merge@meta.data$orig.ident == 2)] <- 14
age[which(merge@meta.data$orig.ident == 3)] <- 19
age[which(merge@meta.data$orig.ident == 4)] <- 13
age[which(merge@meta.data$orig.ident == 5)] <- 19
age[which(merge@meta.data$orig.ident == 6)] <- 22
age[which(merge@meta.data$orig.ident == 7)] <- 4
age[which(merge@meta.data$orig.ident == 8)] <- 12
age[which(merge@meta.data$orig.ident == 9)] <- 6
age[which(merge@meta.data$orig.ident == 10)] <- 21
age[which(merge@meta.data$orig.ident == 11)] <- 19
age[which(merge@meta.data$orig.ident == 12)] <- 44
age[which(merge@meta.data$orig.ident == 13)] <- 34
age[which(merge@meta.data$orig.ident == 14)] <- 21
age[which(merge@meta.data$orig.ident == 15)] <- 19
age[which(merge@meta.data$orig.ident == 16)] <- 46
age[which(merge@meta.data$orig.ident == 17)] <- 24
age[which(merge@meta.data$orig.ident == 18)] <- 54
age[which(merge@meta.data$orig.ident == 19)] <- 39
age[which(merge@meta.data$orig.ident == 20)] <- 33
age[which(merge@meta.data$orig.ident == 21)] <- 27
age[which(merge@meta.data$orig.ident == 22)] <- 25
age[which(merge@meta.data$orig.ident == 23)] <- 49
AddMetaData(merge,metadata=age,col.name="age") -> merge

sex <- rep("male",ncol(merge))

sex[which(merge@meta.data$orig.ident == 2)] <- "female"
sex[which(merge@meta.data$orig.ident == 5)] <- "female"
sex[which(merge@meta.data$orig.ident == 7)] <- "female"
sex[which(merge@meta.data$orig.ident == 12)] <- "female"
sex[which(merge@meta.data$orig.ident == 13)] <- "female"
sex[which(merge@meta.data$orig.ident == 14)] <- "female"
sex[which(merge@meta.data$orig.ident == 16)] <- "female"
sex[which(merge@meta.data$orig.ident == 18)] <- "female"
sex[which(merge@meta.data$orig.ident == 21)] <- "female"
AddMetaData(merge,metadata=sex,col.name="sex") -> merge


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


