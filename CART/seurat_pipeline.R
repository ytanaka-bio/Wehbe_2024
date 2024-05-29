library(Seurat)
library(Matrix)

#CSP <- Read10X_h5("GSM5502393_CTRL_0h_CSP_raw_feature_bc_matrix.h5")
#RNA <- Read10X_h5("GSM5502394_CTRL_0h_RNA_raw_feature_bc_matrix.h5")
#colnames(RNA) <- paste(colnames(RNA),"_Ctrl0h",sep="")
#colnames(CSP) <- paste(colnames(CSP),"_Ctrl0h",sep="")

#merge <- rbind(RNA,CSP)

CSP <- Read10X_h5("GSM5502395_CTRL_24h_CSP_raw_feature_bc_matrix.h5")
RNA <- Read10X_h5("GSM5502396_CTRL_24h_RNA_raw_feature_bc_matrix.h5")
colnames(RNA) <- paste(colnames(RNA),"_Ctrl24h",sep="")
colnames(CSP) <- paste(colnames(CSP),"_Ctrl24h",sep="")

merge <- rbind(RNA,CSP)
protein <- rownames(CSP)

CSP <- Read10X_h5("GSM5502397_STIM_24h_CSP_raw_feature_bc_matrix.h5")
RNA <- Read10X_h5("GSM5502398_STIM_24h_RNA_raw_feature_bc_matrix.h5")
colnames(RNA) <- paste(colnames(RNA),"_STIM24h",sep="")
colnames(CSP) <- paste(colnames(CSP),"_STIM24h",sep="")

merge <- cbind(merge,rbind(RNA,CSP))
data <- merge

merge <- CreateSeuratObject(counts=data,project="combined",names.field = 2, names.delim = "_")
mito.gene <-  rownames(merge)[grep("MT-",rownames(merge))]
percent.mito <- Matrix::colSums(merge@assays$RNA[mito.gene, ])/Matrix::colSums(merge@assays$RNA)
AddMetaData(merge,metadata=percent.mito,col.name="percent.mito") -> merge

treat <- rep("control",ncol(merge))
treat[which(merge@meta.data$orig.ident == "STIM24h")] <- "treat"
AddMetaData(merge,metadata=treat,col.name="Treat") -> merge

#merge <- subset(merge,subset= nFeature_RNA > 100 & nFeature_RNA < 7000 & nCount_RNA > 500 & nCount_RNA < 20000 & percent.mito < 0.05)

#RNA and Protein separately
cells <- colnames(merge)
colSums(data[grep("TotalSeqC",rownames(data))[1:3],cells]) -> HTO_count
colSums(data[grep("TotalSeqC",rownames(data))[-1:-3],cells]) -> ADT_count
colSums(data[rownames(RNA)[-1*grep("MT-",rownames(RNA))],cells]) -> RNA_count
cells <- cells[which(HTO_count != 0 & ADT_count != 0 & RNA_count!=0)]
#RNA_count <- RNA_count[cells]
#ADT_count <- ADT_count[cells]
#HTO_count <- HTO_count[cells]
#background <- cells[which(log10(RNA_count) < 2.75 & log10(ADT_count) > 2 & log10(ADT_count) < 2.9)]


merge_RNA <- CreateSeuratObject(counts=data[rownames(RNA),cells],project="combined",names.field = 2, names.delim = "_")
merge_CSP <- CreateSeuratObject(counts=data[grep("TotalSeqC",rownames(data))[-1:-3],cells],project="combined",names.field = 2, names.delim = "_")

merge_HTO <- CreateSeuratObject(counts=data[grep("TotalSeqC",rownames(data))[1:3],cells],project="combined",names.field = 2, names.delim = "_")


####### not used ##########
#merge
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

source("~/scratch/snorkel_v4/weasel/github/sgMatrix_table.R")
sgMatrix_table(merge.integrated,"matrix/")

############################

#merge protein
library(dsb)

merge_CSP <- NormalizeData(merge_CSP)
merge_CSP <- FindVariableFeatures(merge_CSP,selection.method = "mean.var.plot")
merge_CSP <- ScaleData(merge_CSP,features=VariableFeatures(merge_CSP))
merge_CSP[['HTO']] <- CreateAssayObject(counts=merge_HTO@assays$RNA@counts)

merge_CSP <- NormalizeData(merge_CSP,assay="HTO",normalization.method="CLR")
merge_CSP <- HTODemux(merge_CSP,assay="HTO",positive.quantile=0.99)

singlet <- colnames(merge_CSP)[which(merge_CSP@meta.data$HTO_classification.global == "Singlet")]
negative <- colnames(merge_CSP)[which(merge_CSP@meta.data$HTO_classification.global == "Negative")]

isotype.controls <- rownames(merge_CSP)[29:31]
cells.dsb.norm = DSBNormalizeProtein(cell_protein_matrix=merge_CSP@assays$RNA@counts[,singlet],empty_drop_matrix=merge_CSP@assays$RNA@counts[,negative],denoise.counts=T,use.isotype.control=T,isotype.control.name.vec=isotype.controls)

merge_RNA <- subset(merge_RNA,cells=singlet)


#merge_RNA
merge_RNA.list <- SplitObject(object=merge_RNA,split.by="orig.ident")

for(i in 1:length(merge_RNA.list)){
merge_RNA.list[[i]] <- NormalizeData(merge_RNA.list[[i]],verbose=F)
merge_RNA.list[[i]] <- FindVariableFeatures(merge_RNA.list[[i]],selection.method="vst",nfeatures=2000,verbose=F)
}

merge_RNA.anchors <- FindIntegrationAnchors(merge_RNA.list,dims=1:20)
merge_RNA.integrated <- IntegrateData(anchorset = merge_RNA.anchors,dims=1:20)

library(ggplot2)
library(cowplot)

DefaultAssay(object = merge_RNA.integrated) <- "integrated"
merge_RNA.integrated <- ScaleData(object = merge_RNA.integrated, verbose = FALSE)
merge_RNA.integrated <- RunPCA(object = merge_RNA.integrated, npcs = 20, verbose = FALSE)
merge_RNA.integrated <- RunUMAP(object = merge_RNA.integrated, reduction = "pca", dims = 1:20)

merge_RNA.integrated <- FindNeighbors(merge_RNA.integrated,dims=1:20,reduction="pca")
merge_RNA.integrated <- FindClusters(merge_RNA.integrated)
merge_RNA.integrated@meta.data$Treat <- merge.integrated@meta.data[colnames(merge_RNA.integrated),"Treat"]
merge_RNA.integrated[['ADT']] <- CreateAssayObject(data = cells.dsb.norm)

save(merge_RNA.integrated,file="merge_RNA.dat")
save(merge_RNA.list,file="merge_RNA_list.dat")
save(merge_RNA.anchors,file="merge_RNA_anchors.dat")

source("/home/ytanaka/scratch/snorkel_v4/weasel/github/sgMatrix_table_CITE.R")
sgMatrix_table_CITE(merge_RNA.integrated,"matrix3")

#merge RNA CART
merge_RNA_CART <- subset(merge_RNA,cells=c(cart_bc_ctrl,cart_bc_stim))
merge_RNA_CART.list <- SplitObject(object=merge_RNA_CART,split.by="orig.ident")

for(i in 1:length(merge_RNA_CART.list)){
merge_RNA_CART.list[[i]] <- NormalizeData(merge_RNA_CART.list[[i]],verbose=F)
merge_RNA_CART.list[[i]] <- FindVariableFeatures(merge_RNA_CART.list[[i]],selection.method="vst",nfeatures=2000,verbose=F)
}

merge_RNA_CART.anchors <- FindIntegrationAnchors(merge_RNA_CART.list,dims=1:20)
merge_RNA_CART.integrated <- IntegrateData(anchorset = merge_RNA_CART.anchors,dims=1:20)

library(ggplot2)
library(cowplot)

DefaultAssay(object = merge_RNA_CART.integrated) <- "integrated"
merge_RNA_CART.integrated <- ScaleData(object = merge_RNA_CART.integrated, verbose = FALSE)
merge_RNA_CART.integrated <- RunPCA(object = merge_RNA_CART.integrated, npcs = 20, verbose = FALSE)
merge_RNA_CART.integrated <- RunUMAP(object = merge_RNA_CART.integrated, reduction = "pca", dims = 1:20)

merge_RNA_CART.integrated <- FindNeighbors(merge_RNA_CART.integrated,dims=1:20,reduction="pca")
merge_RNA_CART.integrated <- FindClusters(merge_RNA_CART.integrated)
merge_RNA_CART.integrated@meta.data$Treat <- merge.integrated@meta.data[colnames(merge_RNA_CART.integrated),"Treat"]
merge_RNA_CART.integrated[['ADT']] <- CreateAssayObject(data = cells.dsb.norm[1:28,c(cart_bc_ctrl,cart_bc_stim)])

save(merge_RNA_CART.integrated,file="merge_RNA_CART.dat")
save(merge_RNA_CART.list,file="merge_RNA_CART_list.dat")
save(merge_RNA_CART.anchors,file="merge_RNA_CART_anchors.dat")

sgMatrix_table_CITE(merge_RNA_CART.integrated,"matrix4")

#CD4 and CD8
merge_RNA_CART_CD4 <- subset(merge_RNA_CART, cells = colnames(merge_RNA_CART_CD4.integrated))
merge_RNA_CART_CD8 <- subset(merge_RNA_CART, cells = colnames(merge_RNA_CART_CD8.integrated))

merge_RNA_CART_CD4.list <- SplitObject(object=merge_RNA_CART_CD4,split.by="orig.ident")

for(i in 1:length(merge_RNA_CART_CD4.list)){
merge_RNA_CART_CD4.list[[i]] <- NormalizeData(merge_RNA_CART_CD4.list[[i]],verbose=F)
merge_RNA_CART_CD4.list[[i]] <- FindVariableFeatures(merge_RNA_CART_CD4.list[[i]],selection.method="vst",nfeatures=2000,verbose=F)
}

merge_RNA_CART_CD4.anchors <- FindIntegrationAnchors(merge_RNA_CART_CD4.list,dims=1:20)
merge_RNA_CART_CD4.integrated <- IntegrateData(anchorset = merge_RNA_CART_CD4.anchors,dims=1:20)

library(ggplot2)
library(cowplot)

DefaultAssay(object = merge_RNA_CART_CD4.integrated) <- "integrated"
merge_RNA_CART_CD4.integrated <- ScaleData(object = merge_RNA_CART_CD4.integrated, verbose = FALSE)
merge_RNA_CART_CD4.integrated <- RunPCA(object = merge_RNA_CART_CD4.integrated, npcs = 20, verbose = FALSE)
merge_RNA_CART_CD4.integrated <- RunUMAP(object = merge_RNA_CART_CD4.integrated, reduction = "pca", dims = 1:20)

merge_RNA_CART_CD4.integrated <- FindNeighbors(merge_RNA_CART_CD4.integrated,dims=1:20,reduction="pca")
merge_RNA_CART_CD4.integrated <- FindClusters(merge_RNA_CART_CD4.integrated)
merge_RNA_CART_CD4.integrated@meta.data$Treat <- merge.integrated@meta.data[colnames(merge_RNA_CART_CD4.integrated),"Treat"]
merge_RNA_CART_CD4.integrated[['ADT']] <- CreateAssayObject(data = cells.dsb.norm[1:28,colnames(merge_RNA_CART_CD4.integrated)])

save(merge_RNA_CART_CD4.integrated,file="merge_RNA_CART_CD4.dat")
save(merge_RNA_CART_CD4.list,file="merge_RNA_CART_CD4_list.dat")
save(merge_RNA_CART_CD4.anchors,file="merge_RNA_CART_CD4_anchors.dat")


merge_RNA_CART_CD8.list <- SplitObject(object=merge_RNA_CART_CD8,split.by="orig.ident")

for(i in 1:length(merge_RNA_CART_CD8.list)){
merge_RNA_CART_CD8.list[[i]] <- NormalizeData(merge_RNA_CART_CD8.list[[i]],verbose=F)
merge_RNA_CART_CD8.list[[i]] <- FindVariableFeatures(merge_RNA_CART_CD8.list[[i]],selection.method="vst",nfeatures=2000,verbose=F)
}

merge_RNA_CART_CD8.anchors <- FindIntegrationAnchors(merge_RNA_CART_CD8.list,dims=1:20)
merge_RNA_CART_CD8.integrated <- IntegrateData(anchorset = merge_RNA_CART_CD8.anchors,dims=1:20)

library(ggplot2)
library(cowplot)

DefaultAssay(object = merge_RNA_CART_CD8.integrated) <- "integrated"
merge_RNA_CART_CD8.integrated <- ScaleData(object = merge_RNA_CART_CD8.integrated, verbose = FALSE)
merge_RNA_CART_CD8.integrated <- RunPCA(object = merge_RNA_CART_CD8.integrated, npcs = 20, verbose = FALSE)
merge_RNA_CART_CD8.integrated <- RunUMAP(object = merge_RNA_CART_CD8.integrated, reduction = "pca", dims = 1:20)

merge_RNA_CART_CD8.integrated <- FindNeighbors(merge_RNA_CART_CD8.integrated,dims=1:20,reduction="pca")
merge_RNA_CART_CD8.integrated <- FindClusters(merge_RNA_CART_CD8.integrated)
merge_RNA_CART_CD8.integrated@meta.data$Treat <- merge.integrated@meta.data[colnames(merge_RNA_CART_CD8.integrated),"Treat"]
merge_RNA_CART_CD8.integrated[['ADT']] <- CreateAssayObject(data = cells.dsb.norm[1:28,colnames(merge_RNA_CART_CD8.integrated)])

save(merge_RNA_CART_CD8.integrated,file="merge_RNA_CART_CD8.dat")
save(merge_RNA_CART_CD8.list,file="merge_RNA_CART_CD8_list.dat")
save(merge_RNA_CART_CD8.anchors,file="merge_RNA_CART_CD8_anchors.dat")


#CD4 and CD8
df <- data.frame(CD4=cells.dsb.norm[1,colnames(merge_RNA)],CD8=cells.dsb.norm[2,colnames(merge_RNA)])
df[df < -2] <- -2
df[df[,1] > 8,1] <- 8
df[df[,2] > 7,2] <- 7
library(fields)
library(gplots)
h2d <- hist2d(df[,1],df[,2],nbin=100,show=F)
samp <- log10(h2d$counts+1)
samp[samp > 1.5] <- 1.5
image.plot(h2d$x,h2d$y,samp,col=c("white",tim.colors(127)))
samp <- cells.dsb.norm[1:2,colnames(merge_RNA)]

merge_RNA_CD4 <- subset(merge_RNA,cells=CD4)
merge_RNA_CD8 <- subset(merge_RNA,cells=CD8)
merge_RNA_DN <- subset(merge_RNA,cells=DN)

#merge_RNA_CD4
merge_RNA_CD4.list <- SplitObject(object=merge_RNA_CD4,split.by="orig.ident")

for(i in 1:length(merge_RNA_CD4.list)){
merge_RNA_CD4.list[[i]] <- NormalizeData(merge_RNA_CD4.list[[i]],verbose=F)
merge_RNA_CD4.list[[i]] <- FindVariableFeatures(merge_RNA_CD4.list[[i]],selection.method="vst",nfeatures=2000,verbose=F)
}

merge_RNA_CD4.anchors <- FindIntegrationAnchors(merge_RNA_CD4.list,dims=1:20)
merge_RNA_CD4.integrated <- IntegrateData(anchorset = merge_RNA_CD4.anchors,dims=1:20)

library(ggplot2)
library(cowplot)

DefaultAssay(object = merge_RNA_CD4.integrated) <- "integrated"
merge_RNA_CD4.integrated <- ScaleData(object = merge_RNA_CD4.integrated, verbose = FALSE)
merge_RNA_CD4.integrated <- RunPCA(object = merge_RNA_CD4.integrated, npcs = 20, verbose = FALSE)
merge_RNA_CD4.integrated <- RunUMAP(object = merge_RNA_CD4.integrated, reduction = "pca", dims = 1:20)

merge_RNA_CD4.integrated <- FindNeighbors(merge_RNA_CD4.integrated,dims=1:20,reduction="pca")
merge_RNA_CD4.integrated <- FindClusters(merge_RNA_CD4.integrated)
merge_RNA_CD4.integrated@meta.data$Treat <- merge.integrated@meta.data[colnames(merge_RNA_CD4.integrated),"Treat"]
merge_RNA_CD4.integrated[['ADT']] <- CreateAssayObject(data = cells.dsb.norm[,colnames(merge_RNA_CD4.integrated)])

save(merge_RNA_CD4.integrated,file="merge_RNA_CD4.dat")
save(merge_RNA_CD4.list,file="merge_RNA_CD4_list.dat")
save(merge_RNA_CD4.anchors,file="merge_RNA_CD4_anchors.dat")

#merge_RNA_CD8
merge_RNA_CD8.list <- SplitObject(object=merge_RNA_CD8,split.by="orig.ident")

for(i in 1:length(merge_RNA_CD8.list)){
merge_RNA_CD8.list[[i]] <- NormalizeData(merge_RNA_CD8.list[[i]],verbose=F)
merge_RNA_CD8.list[[i]] <- FindVariableFeatures(merge_RNA_CD8.list[[i]],selection.method="vst",nfeatures=2000,verbose=F)
}

merge_RNA_CD8.anchors <- FindIntegrationAnchors(merge_RNA_CD8.list,dims=1:20)
merge_RNA_CD8.integrated <- IntegrateData(anchorset = merge_RNA_CD8.anchors,dims=1:20)

library(ggplot2)
library(cowplot)

DefaultAssay(object = merge_RNA_CD8.integrated) <- "integrated"
merge_RNA_CD8.integrated <- ScaleData(object = merge_RNA_CD8.integrated, verbose = FALSE)
merge_RNA_CD8.integrated <- RunPCA(object = merge_RNA_CD8.integrated, npcs = 20, verbose = FALSE)
merge_RNA_CD8.integrated <- RunUMAP(object = merge_RNA_CD8.integrated, reduction = "pca", dims = 1:20)

merge_RNA_CD8.integrated <- FindNeighbors(merge_RNA_CD8.integrated,dims=1:20,reduction="pca")
merge_RNA_CD8.integrated <- FindClusters(merge_RNA_CD8.integrated)
merge_RNA_CD8.integrated@meta.data$Treat <- merge.integrated@meta.data[colnames(merge_RNA_CD8.integrated),"Treat"]
merge_RNA_CD8.integrated[['ADT']] <- CreateAssayObject(data = cells.dsb.norm[,colnames(merge_RNA_CD8.integrated)])

save(merge_RNA_CD8.integrated,file="merge_RNA_CD8.dat")
save(merge_RNA_CD8.list,file="merge_RNA_CD8_list.dat")
save(merge_RNA_CD8.anchors,file="merge_RNA_CD8_anchors.dat")

#merge_RNA_DN
merge_RNA_DN.list <- SplitObject(object=merge_RNA_DN,split.by="orig.ident")

for(i in 1:length(merge_RNA_DN.list)){
merge_RNA_DN.list[[i]] <- NormalizeData(merge_RNA_DN.list[[i]],verbose=F)
merge_RNA_DN.list[[i]] <- FindVariableFeatures(merge_RNA_DN.list[[i]],selection.method="vst",nfeatures=2000,verbose=F)
}

merge_RNA_DN.anchors <- FindIntegrationAnchors(merge_RNA_DN.list,dims=1:20)
merge_RNA_DN.integrated <- IntegrateData(anchorset = merge_RNA_DN.anchors,dims=1:20)

library(ggplot2)
library(cowplot)

DefaultAssay(object = merge_RNA_DN.integrated) <- "integrated"
merge_RNA_DN.integrated <- ScaleData(object = merge_RNA_DN.integrated, verbose = FALSE)
merge_RNA_DN.integrated <- RunPCA(object = merge_RNA_DN.integrated, npcs = 20, verbose = FALSE)
merge_RNA_DN.integrated <- RunUMAP(object = merge_RNA_DN.integrated, reduction = "pca", dims = 1:20)

merge_RNA_DN.integrated <- FindNeighbors(merge_RNA_DN.integrated,dims=1:20,reduction="pca")
merge_RNA_DN.integrated <- FindClusters(merge_RNA_DN.integrated)
merge_RNA_DN.integrated@meta.data$Treat <- merge.integrated@meta.data[colnames(merge_RNA_DN.integrated),"Treat"]
merge_RNA_DN.integrated[['ADT']] <- CreateAssayObject(data = cells.dsb.norm[,colnames(merge_RNA_DN.integrated)])

save(merge_RNA_DN.integrated,file="merge_RNA_DN.dat")
save(merge_RNA_DN.list,file="merge_RNA_DN_list.dat")
save(merge_RNA_DN.anchors,file="merge_RNA_DN_anchors.dat")

sgMatrix_table_CITE(merge_RNA_CD4.integrated,"matrix_CD4_all")
sgMatrix_table_CITE(merge_RNA_CD8.integrated,"matrix_CD8_all")
sgMatrix_table_CITE(merge_RNA_DN.integrated,"matrix_DN_all")
