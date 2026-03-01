library(Seurat)
library(UCell)
library(clusterProfiler)
library(progressr)
library(Matrix)
handlers(global = TRUE)

setwd("/slurm/home/yrd/liaolab/nieshuyang/AHBA_sc")
seurat_obj <- readRDS("/scRNA/adult_ctx.rds")
# Load gene specificity score matrix in R
gss <- readMM("/scRNA/gss.mtx")
seurat_obj@assays[["RNA"]]@layers[["gss"]] <- t(gss)
rownames(seurat_obj@assays[["RNA"]]@layers[["data"]]) <- rownames(seurat_obj)
colnames(seurat_obj@assays[["RNA"]]@layers[["data"]]) <- colnames(seurat_obj)
rownames(seurat_obj@assays[["RNA"]]@layers[["counts"]]) <- rownames(seurat_obj)
colnames(seurat_obj@assays[["RNA"]]@layers[["counts"]]) <- colnames(seurat_obj)
rownames(seurat_obj@assays[["RNA"]]@layers[["gss"]]) <- rownames(seurat_obj)
colnames(seurat_obj@assays[["RNA"]]@layers[["gss"]]) <- colnames(seurat_obj)
# Load gene set
geneset <- read.gmt("/GMT/GO_Biological_Process_2025.gmt")
geneset <- split(geneset$gene, geneset$term)
geneset <- lapply(geneset, function(x) intersect(x, rownames(seurat_obj)))
pathways <- c(
  "Cellular Respiration (GO:0045333)",
  "Oxidative Phosphorylation (GO:0006119)",
  "Tricarboxylic Acid Cycle (GO:0006099)",
  "Mitochondrial ATP Synthesis Coupled Electron Transport (GO:0042775)",
  "Energy Derivation by Oxidation of Organic Compounds (GO:0015980)",
  "Fatty Acid Beta-Oxidation (GO:0006635)",
  "Central Nervous System Myelination (GO:0022010)"
)
geneset_subset <- geneset[names(geneset) %in% pathways]
# Score cells through UCell
with_progress({ 
  seurat_obj <- AddModuleScore_UCell(seurat_obj, features=geneset_subset, slot="gss", name=NULL)
})
ucell_res <- seurat_obj@meta.data
ucell_res <- ucell_res[,colnames(ucell_res) %in% pathways]
write.csv(ucell_res,file="/GMT/UCell_scores.csv")