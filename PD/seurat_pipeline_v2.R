#module load R/3.5.0-foss-2016b-avx2
library(Seurat)
library(Matrix)
library(EnsDb.Hsapiens.v86)
library(GenomicRanges)
library(Signac)

list <- read.table("RNA_h5.txt",sep="\t",header=F)
list <- as.character(list[,1])
list2 <- read.table("frag_list.txt",sep="\t",header=F)
list2 <- as.character(list2[,1])
id <- read.table("ID_list.txt",sep="\t",header=F)
id <- as.character(id[,1])
orig_name <- vector("list",length(list))
names(orig_name) <- list

data_stat <- vector("list",4)
names(data_stat) <- c("nCount_RNA","nFeature_RNA","nCount_ATAC","nFeature_ATAC")
for(i in 1:4){
      data_stat[[i]] <- vector("list",length(list))
      names(data_stat[[i]]) <- id
}
filtered_cells <- matrix(0,length(id),2)
rownames(filtered_cells) <- id
colnames(filtered_cells) <- c("Post","Pre")

#QC
for(i in 1:length(list)){
all <- Read10X_h5(list[i])
data <- CreateSeuratObject(counts=all[[1]])
data[['ATAC']] <- CreateChromatinAssay(counts=all[[2]],sep=c(":","-"),genome="hg38",fragments=list2[i])

data_stat[[1]][[i]] <- data@meta.data$nCount_RNA
data_stat[[2]][[i]] <- data@meta.data$nFeature_RNA
data_stat[[3]][[i]] <- data@meta.data$nCount_ATAC
data_stat[[4]][[i]] <- data@meta.data$nFeature_ATAC
#}
#save(data_stat,file="data_stat.dat")
filtered_cells[i,2] <- ncol(data)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
Annotation(data[["ATAC"]]) <- annotations

#QC-ATAC
DefaultAssay(data) <- "ATAC"
data <- TSSEnrichment(data)
data <- NucleosomeSignal(data)
data$blacklist_fraction <- FractionCountsInRegion(object = data,assay = 'ATAC',regions = blacklist_hg38)
data <- subset(x = data,subset = blacklist_fraction < 0.02 & TSS.enrichment > 2 & nucleosome_signal < 4 & nCount_ATAC > 1000 & nCount_ATAC < 60000 & nFeature_ATAC > 1000 & nFeature_ATAC < 25000)

#QC RNA
DefaultAssay(data) <- "RNA"
mito.gene <-  rownames(data)[grep("MT-",rownames(data))]
percent.mito <- Matrix::colSums(data@assays$RNA[mito.gene, ])/Matrix::colSums(data@assays$RNA)
AddMetaData(data,metadata=percent.mito,col.name="percent.mito") -> data
data <- subset(data,subset= nFeature_RNA > 100 & nFeature_RNA < 7000 & nCount_RNA > 500 & nCount_RNA < 20000 & percent.mito < 0.05)
orig_name[[i]] <- colnames(data)
filtered_cells[i,1] <- ncol(data)
}
save(data_stat,file="data_stat.dat")
save(orig_name,file="orig_name_v2.dat")
save(filtered_cells,file="filtered_cells.dat")

#library filter
use_lib <- which(filtered_cells[,1]/filtered_cells[,2] > 0.4)
list[use_lib] -> list
list2[use_lib] -> list2
id[use_lib] -> id
orig_name2 <- vector("list",length(list))
names(orig_name2) <- list

for(i in 1:length(list)){
      orig_name2[[list[i]]] <- orig_name[[list[i]]]
}
orig_name <- orig_name2
rm(orig_name2)
save(orig_name,file="orig_name_v2_2.dat")
save(use_lib,file="use_lib_v2.dat")

#Normalize
all <- Read10X_h5(list[1])
all_peak <- vector("list",length(list)) 
names(all_peak) <- list
all[[1]] <- all[[1]][,orig_name[[1]]]
all[[2]] <- all[[2]][,orig_name[[1]]]
colnames(all[[1]]) <- paste(colnames(all[[1]]),use_lib[1],sep="_")
colnames(all[[2]]) <- paste(colnames(all[[2]]),use_lib[1],sep="_")
all_peak[[1]] <- all[[2]]
all <- all[[1]]

for(i in 2:length(list)){
      data <- Read10X_h5(list[i])
      data[[1]] <- data[[1]][,orig_name[[i]]]
      data[[2]] <- data[[2]][,orig_name[[i]]]
      colnames(data[[1]]) <- paste(colnames(data[[1]]),use_lib[i],sep="_")
      colnames(data[[2]]) <- paste(colnames(data[[2]]),use_lib[i],sep="_")
      all <- cbind(all,data[[1]])
      all_peak[[i]] <- data[[2]]
}
rm(data)

merge <- CreateSeuratObject(counts=all,project="combined",names.field = 2, names.delim = "_")

mito.gene <-  rownames(merge)[grep("MT-",rownames(merge))]
percent.mito <- Matrix::colSums(merge@assays$RNA[mito.gene, ])/Matrix::colSums(merge@assays$RNA)
AddMetaData(merge,metadata=percent.mito,col.name="percent.mito") -> merge
rm(all)


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
merge.integrated <- RunPCA(object = merge.integrated, npcs = 30, verbose = FALSE)
merge.integrated <- RunUMAP(object = merge.integrated, reduction = "pca", dims = 1:30)

merge.integrated <- FindNeighbors(merge.integrated,dims=1:20,reduction="pca")
merge.integrated <- FindClusters(merge.integrated)

save(merge.integrated,file="merge_v2.dat")
save(merge.list,file="merge_list_v2.dat")
save(merge.anchors,file="merge_anchors_v2.dat")


#ATAC combine peak
library(Signac)
library(GenomicRanges)
list <- read.table("ATAC_bed.txt",sep="\t",header=F)
list <- as.character(list[,1])[use_lib]
peak_vec <- vector("list",length(list))
names(peak_vec) <- list
for(i in 1:length(list)){
peak_vec[[i]] <- read.table(list[i],sep="\t",col.names=c("chr", "start", "end"))
}

peak_vec -> gr_vec
gr_list <- makeGRangesFromDataFrame(peak_vec[[1]])
for(i in 2:length(list)){
gr_vec[[i]] <- makeGRangesFromDataFrame(peak_vec[[i]])
gr_list <- c(gr_list,gr_vec[[i]])
}
#GRangesList(gr_vec) -> gr_list
combined.peaks <- reduce(x=gr_list)
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
save(combined.peaks,file="combined.peaks_v2.dat")

#ATAC fragment
list <- read.table("frag_list.txt",sep="\t",header=F)
list <- as.character(list[,1])[use_lib]
frag_vec <- vector("list",length(list))
names(frag_vec) <- list
for(i in 1:length(list)){
frag_vec[[i]] <- CreateFragmentObject(path=list[i],cells=orig_name[[i]])
}
save(frag_vec,file="frag_vec_v2.dat")

#Feature matrix
f_mat <- frag_vec
for(i in 1:length(list)){
f_mat[[i]] <- FeatureMatrix(fragments=frag_vec[[i]],features=combined.peaks,cells=orig_name[[i]])
}
for(i in 1:length(list)){
colnames(f_mat[[i]]) <- paste(colnames(f_mat[[i]]),use_lib[i],sep="_")
}
save(f_mat,file="f_mat_v2.dat")

#combine all
all_mat <- f_mat[[1]]
for(i in 2:length(f_mat)){
all_mat <- cbind(all_mat,f_mat[[i]])
}
C_obj <- CreateChromatinAssay(all_mat,genome="hg38",assay="ATAC")
C_obj <- CreateSeuratObject(C_obj,assay="ATAC",names.field = 2, names.delim="_")
save(C_obj,file="C_obj_combined_v2.dat")

#Merge peaks
library(harmony)
C_obj <- RunTFIDF(C_obj)
C_obj <- FindTopFeatures(C_obj, min.cutoff = 50)
C_obj <- RunSVD(C_obj)
merge.peaks <- RunHarmony(C_obj,"orig.ident",reduction="lsi",assay.use="ATAC",project.dim=F)
merge.peaks <- RunUMAP(merge.peaks, dims = 2:20, reduction = 'harmony')

#add fragment object
list <- read.table("new_frag_list.txt",sep="\t",header=F)
list <- as.character(list[,1])[use_lib]
new_frag_vec <- vector("list",length(list))
names(new_frag_vec) <- list
for(i in 1:length(list)){
new_frag_vec[[i]] <- CreateFragmentObject(path=list[i],cells=colnames(f_mat[[i]]))
}
save(new_frag_vec,file="new_frag_vec_v2.dat")
Fragments(merge.peaks) <- new_frag_vec
save(merge.peaks,file="merge.peaks_v2.dat")

#Create Object
#C_obj <- frag_vec
#for(i in 1:length(list)){
#C_obj[[i]] <- CreateChromatinAssay(counts = f_mat[[i]], fragments = frag_vec[[i]], genome="hg38",assay="ATAC")
#C_obj[[i]] <- RunTFIDF(C_obj[[i]])
#C_obj[[i]] <- CreateSeuratObject(C_obj[[i]], assay = "ATAC",names.field = 2, names.delim = "_")
#C_obj[[i]]$orig.ident <- as.character(i)
#}
#save(C_obj,file="C_obj_v2.dat")

#combine object
#peak.anchors <- FindIntegrationAnchors(object.list = C_obj, anchor.features = rownames(C_obj[[1]]),dims=1:20,assay=rep("ATAC",length(C_obj)),k.filter=NA)
#merge.peaks <- IntegrateData(anchorset = peak.anchors,dims = 2:20,preserve.order = TRUE)
#merge.peaks <- RunSVD(merge.peaks,reduction.name = 'integratedLSI" )
#merge.peaks <- RunUMAP(merge.peaks, dims = 2:20, reduction = 'integratedLSI')
#save(merge.peaks,file="merge.peaks_v2.dat")
#save(peak.anchors,file="peak.anchors_v2.dat")

