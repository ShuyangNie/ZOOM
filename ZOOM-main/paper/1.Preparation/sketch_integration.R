## Preprocess and integrate scRNA-seq data ###
# Load packages
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggrepel)
# Set random seed for reproducibility
set.seed(123)
# Avoid memory leak
options(future.globals.maxSize = 1024*3e+09)

## Basic Seurat pre-process pipeline

setwd('/slurm/home/yrd/liaolab/nieshuyang/AHBA_sc')
# Read SeuratObject, or you can construct the SeuratObject de novo
seurat_obj <- readRDS('scRNA/adult_ctx.rds')

# Quality control
# Find MT genes
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
# The threshold was determined by the distribution
seurat_obj <- subset(
    seurat_obj, 
    subset = nFeature_RNA > 250 & nFeature_RNA < 2500 & percent.mt < 5
)

# Allign genes in scRNA-seq dataset with AHBA dataset
ahba_genes <- read.csv("HCPMMP/DS_HCPMMP.csv", row.names=1)
ahba_genes <- rownames(ahba_genes)
seurat_obj <- subset(seurat_obj, features=ahba_genes)

# Normalize data
seurat_obj <- NormalizeData(
    seurat_obj, 
    normalization.method = "LogNormalize", 
    scale.factor = 10000
)

# Make sure rownames and colnames present
rownames(seurat_obj@assays[["RNA"]]@layers[["data"]]) <- rownames(seurat_obj)
colnames(seurat_obj@assays[["RNA"]]@layers[["data"]]) <- colnames(seurat_obj)
rownames(seurat_obj@assays[["RNA"]]@layers[["counts"]]) <- rownames(seurat_obj)
colnames(seurat_obj@assays[["RNA"]]@layers[["counts"]]) <- colnames(seurat_obj)


### Sketch integration pipeline
# Split dataset by batch key that present in `seurat_obj@meta.data`, as controlled by argument `f`
seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]], f=seurat_obj$Dataset)
# Find high variable genes
seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
# Sample representative cells from each dataset
seurat_obj <- SketchData(
    seurat_obj, 
    ncell = 5000, # Or set a lower number
    method = "LeverageScore", 
    sketched.assay = "sketch"
)
# Perform integration on the sketched cells across samples
DefaultAssay(seurat_obj) <- "sketch"

# Find high variable genes, scale the data and run PCA on representative cells
seurat_obj <- FindVariableFeatures(seurat_obj, verbose = F)
seurat_obj <- ScaleData(seurat_obj, verbose = F)
seurat_obj <- RunPCA(seurat_obj, verbose = F)
# Integrate these representative cells
seurat_obj <- IntegrateLayers(
    seurat_obj, 
    method = RPCAIntegration, # You can try other integration methods
    orig = "pca", 
    new.reduction = "integrated.rpca",
    dims = 1:30, 
    k.anchor = 20, # You can try other values
    verbose = F
)

# Cluster the integrated data (representative cells)
seurat_obj <- FindNeighbors(seurat_obj, reduction = "integrated.rpca", dims = 1:30)
seurat_obj <- RunUMAP(
    seurat_obj, 
    reduction = "integrated.rpca", 
    dims = 1:30, 
    return.model = T, 
    verbose = F
)
# Plot the umap to check if representative cells from different batches have been well mixed
DimPlot(
    seurat_obj, 
    group.by = "Dataset",
    reduction = "umap"
)

# Integrate the full datasets
seurat_obj <- ProjectIntegration(
    seurat_obj, 
    sketched.assay = "sketch", 
    assay = "RNA", 
    reduction = "integrated.rpca"
)
seurat_obj <- RunUMAP(
    seurat_obj, 
    reduction = "integrated.rpca.full", 
    dims = 1:30, 
    reduction.name = "umap.full"
)
# Plot the umap to check if representative cells from different batches have been well mixed
DimPlot(
    seurat_obj, 
    reduction = "umap.full", 
    group.by = "Dataset"
)
# Restore the dataset
seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])
DefaultAssay(seurat_obj) <- "RNA"

# Save results as .rds file
saveRDS(seurat_obj, file='scRNA/adult_ctx.rds')