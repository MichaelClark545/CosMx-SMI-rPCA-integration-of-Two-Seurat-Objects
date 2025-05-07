# CosMx-SMI-rPCA-integration-of-Two-Seurat-Objects
Integrate with two Seurat objects


# Set allowed size to 64 gig
options(future.globals.maxSize=64000*1024^2)

#Seurat Objects
library(Seurat)
library(ggplot2)
library(patchwork)
Seurat_Object_1 <- readRDS("Seurat_Object_1_File_Path")
Seurat_Object_2 <- readRDS("Seurat_Object_2_File_Path")

#Merge the two Seurat objects
merged_seurat <- merge(Seurat_Object_1, y = Seurat_Object_2, add.cells.ids = c("sample1", "sample2"))

# split the dataset
# Tissue is listed in metadata of each Seurat object
split_seurat_list <- SplitObject(merged_seurat, split.by = "Tissue")

# Preprocess the Seurat objects
seuratObj.list <- lapply(X = split_seurat_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# PCA Organization
# Prior to starting this step, you will need an extensive amount of RAM to run PCA and integrate the data. Seurat objects with more than 50,000 cells will require at least a 64GB. 
features <- SelectIntegrationFeatures(object.list = seuratObj.list)
seuratObj.list <- lapply(X = seuratObj.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

#Integration
immune.anchors <- FindIntegrationAnchors(object.list = seuratObj.list, anchor.features = features, reduction = "rpca")

#Integration Part 2
immune.combined <- IntegrateData(anchorset = immune.anchors)

#Integration Part 3
DefaultAssay(immune.combined) <- "integrated"
saveRDS(immune.combined, file = "./immune.combined.rds")
#UMAP
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "Tissue")
# The group by function will take metadata from the second group of cell types from the top. The other cell types the metadata gives will not work for the UMAP creation. 
p2 <- DimPlot(immune.combined, reduction = "umap", group.by = "Cell_cluster_in_metadata", label = TRUE, repel = TRUE) 
p1+p2
