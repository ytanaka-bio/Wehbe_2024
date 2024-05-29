library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(scales)
library(gplots)
library(Matrix)
library(ggplot2)


cds <- as.cell_data_set(merge_RNA_CART_CD4.integrated)
cds <- cluster_cells(cds,reduction_method="UMAP")
cds <- learn_graph(cds,use_partition = TRUE)
cds <- order_cells(cds,reduction_method="UMAP",root_cells = colnames(merge_RNA_CART_CD4.integrated)[which(merge_RNA_CART_CD4.integrated@active.ident == "1")])
plot_cells(cds,color_cells_by="pseudotime",show_trajectory_graph = F,cell_size = 1)
type <- rep("a",nrow(pData(cds)))
type[which(pData(cds)$Treat == "control")] <- "b"
pData(cds)$Type <- type
plot_cells(cds,color_cells_by="Type",show_trajectory_graph = F, cell_size = 1)
data <- read.table("~/scratch/snorkel_v11/weasel/github/CITE_train_wsl4_result.csv",sep=",",header=T,row.names=1)
pData(cds)$scIDST <- data[rownames(pData(cds)),1]
plot_cells(cds,color_cells_by="scIDST",show_trajectory_graph = F, cell_size = 1) + scale_color_gradient(low="blue",high="red")
pData(cds)$CellCycle <- colMeans(merge_RNA_CART_CD4.integrated@assays$RNA@data[c("TOP2A","MKI67","E2F1"),])[rownames(pData(cds))]
plot_cells(cds,color_cells_by="CellCycle",show_trajectory_graph = F, cell_size = 1) + scale_color_gradient(low="gray",high="red")
save(cds,file="merge_RNA_CART_CD4_cds.dat")

pseudo <- pseudotime(cds)
exp <- merge_RNA_CART_CD4.integrated@assays$ADT@data
samp <- t(apply(exp[-1:-2,names(sort(pseudo))],1,scale))
loess_smooth <- function(x){
x <- data.frame(X=1:length(x),Y=x)
lo <- loess(Y ~ X, x)
return(predict(lo,x[,1]))
}
samp <- t(apply(samp,1,loess_smooth))
samp[samp > 0.5] <- 0.5
samp[samp < -0.5] <- -0.5
heatmap.2(samp,trace="none",scale="none",col=colorpanel(128,"blue","black","yellow"),density.info="none",Colv=F)
cor.test(pseudotime(cds),colMeans(merge_RNA_CART_CD4.integrated@assays$RNA@data[c("TOP2A","MKI67","E2F1"),]))

pseudo <- pData(cds)$scIDST
names(pseudo) <- rownames(pData(cds))
exp <- merge_RNA_CART_CD4.integrated@assays$ADT@data
samp <- t(apply(exp[-1:-2,names(sort(pseudo))],1,scale))
samp <- t(apply(samp,1,loess_smooth))
samp[samp > 0.5] <- 0.5
samp[samp < -0.5] <- -0.5
heatmap.2(samp,trace="none",scale="none",col=colorpanel(128,"blue","black","yellow"),density.info="none",Colv=F)


cds <- as.cell_data_set(merge_RNA_CART_CD8.integrated)
cds <- cluster_cells(cds,reduction_method="UMAP")
cds <- learn_graph(cds,use_partition = TRUE)
cds <- order_cells(cds,reduction_method="UMAP",root_cells = colnames(merge_RNA_CART_CD8.integrated)[which(merge_RNA_CART_CD8.integrated@active.ident == "5")])
plot_cells(cds,color_cells_by="pseudotime",show_trajectory_graph = F,cell_size = 1)
type <- rep("a",nrow(pData(cds)))
type[which(pData(cds)$Treat == "control")] <- "b"
pData(cds)$Type <- type
plot_cells(cds,color_cells_by="Type",show_trajectory_graph = F, cell_size = 1)
data <- read.table("~/scratch/snorkel_v11/weasel/github/CITE_train_wsl4_result.csv",sep=",",header=T,row.names=1)
pData(cds)$scIDST <- data[rownames(pData(cds)),1]
plot_cells(cds,color_cells_by="scIDST",show_trajectory_graph = F, cell_size = 1) + scale_color_gradient(low="blue",high="red")
pData(cds)$CellCycle <- colMeans(merge_RNA_CART_CD8.integrated@assays$RNA@data[c("TOP2A","MKI67","E2F1"),])[rownames(pData(cds))]
plot_cells(cds,color_cells_by="CellCycle",show_trajectory_graph = F, cell_size = 1) + scale_color_gradient(low="gray",high="red")
save(cds,file="merge_RNA_CART_CD8_cds.dat")

pseudo <- pseudotime(cds)
exp <- merge_RNA_CART_CD8.integrated@assays$ADT@data
samp <- t(apply(exp[-1:-2,names(sort(pseudo))],1,scale))
loess_smooth <- function(x){
x <- data.frame(X=1:length(x),Y=x)
lo <- loess(Y ~ X, x)
return(predict(lo,x[,1]))
}
samp <- t(apply(samp,1,loess_smooth))
samp[samp > 0.5] <- 0.5
samp[samp < -0.5] <- -0.5
heatmap.2(samp,trace="none",scale="none",col=colorpanel(128,"blue","black","yellow"),density.info="none",Colv=F)
cor.test(pseudotime(cds),colMeans(merge_RNA_CART_CD8.integrated@assays$RNA@data[c("TOP2A","MKI67","E2F1"),]))

pseudo <- pData(cds)$scIDST
names(pseudo) <- rownames(pData(cds))
exp <- merge_RNA_CART_CD8.integrated@assays$ADT@data
samp <- t(apply(exp[-1:-2,names(sort(pseudo))],1,scale))
samp <- t(apply(samp,1,loess_smooth))
samp[samp > 0.5] <- 0.5
samp[samp < -0.5] <- -0.5
heatmap.2(samp,trace="none",scale="none",col=colorpanel(128,"blue","black","yellow"),density.info="none",Colv=F)
