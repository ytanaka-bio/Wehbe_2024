library(Seurat)
library(Matrix)

merge <- Read10X_h5("merge/outs/count/filtered_feature_bc_matrix.h5")
merge <- CreateSeuratObject(counts=merge,project="combined",names.field = 2, names.delim = "-")
mito.gene <-  rownames(merge)[grep("MT-",rownames(merge))]
percent.mito <- Matrix::colSums(merge@assays$RNA[mito.gene, ])/Matrix::colSums(merge@assays$RNA)
AddMetaData(merge,metadata=percent.mito,col.name="percent.mito") -> merge

type <- rep("control",ncol(merge))
type[which(merge@meta.data$orig.ident == 3)] <- "DHT"
type[which(merge@meta.data$orig.ident == 4)] <- "DHT"
type[which(merge@meta.data$orig.ident == 5)] <- "EST"
type[which(merge@meta.data$orig.ident == 6)] <- "EST"
AddMetaData(merge,metadata=type,col.name="Type") -> merge


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

