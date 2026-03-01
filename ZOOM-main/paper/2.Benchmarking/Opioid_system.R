## Opioid system decoding: VAM and UCell
library(Seurat)
library(VAM)
library(progress)
library(UCell)

setwd("/slurm/home/yrd/liaolab/nieshuyang/AHBA_sc")
seurat_obj <- readRDS("scRNA/adult_ctx.rds")
rownames(seurat_obj@assays[["RNA"]]@layers[["data"]]) <- rownames(seurat_obj)
colnames(seurat_obj@assays[["RNA"]]@layers[["data"]]) <- colnames(seurat_obj)
rownames(seurat_obj@assays[["RNA"]]@layers[["counts"]]) <- rownames(seurat_obj)
colnames(seurat_obj@assays[["RNA"]]@layers[["counts"]]) <- colnames(seurat_obj)
gc()

# This is a copy of VAM function, we excluded some potential error
VAM_Seurat <- function (seurat.data, gene.set.collection, 
                        center = FALSE, gamma = TRUE) {
  # Extract normalized data expression
  normalized.counts = seurat.data@assays$RNA@layers$data
  colnames(normalized.counts) = colnames(seurat.data)
  
  # Estimate technical variation and assign binary weights
  tech.var.prop = getTechVarPropV5(seurat.data)
  names(tech.var.prop) = rownames(seurat.data)
  gene.weights = rep(1, nrow(normalized.counts))
  names(gene.weights) = rownames(seurat.data)
  
  # Transpose
  gene.expr = Matrix::t(normalized.counts)
  rm(normalized.counts)
  p = ncol(gene.expr)
  n = nrow(gene.expr)
  
  # Extract information
  cell.ids = rownames(gene.expr)
  num.sets = length(gene.set.collection)
  set.names = names(gene.set.collection)
  gene.ids = colnames(gene.expr)
  set.sizes = unlist(lapply(gene.set.collection, length))
  min.set.size = min(set.sizes)
  median.set.size = median(set.sizes)
  results = list()
  results$distance.sq = matrix(0, nrow = n, ncol = num.sets, 
                               dimnames = list(cell.ids, set.names))
  results$cdf.value = matrix(0, nrow = n, ncol = num.sets, 
                             dimnames = list(cell.ids, set.names))
  message("Computing VAM distances for ", num.sets, " gene sets, ", 
          n, " cells and ", p, " genes.")
  message("Min set size: ", min.set.size, ", median size: ", 
          median.set.size)
  
  pb <- progress_bar$new(total = num.sets, format = "  [:bar] :percent :elapsed")
  for (i in 1:num.sets) {
    set.members = gene.set.collection[[i]]
    set.size = set.sizes[i]
    set.exprs = gene.expr[, set.members]
    set.weights = gene.weights[set.members]
    names(tech.var.prop) = colnames(gene.expr)
    vam.results = vam(gene.expr = set.exprs, tech.var.prop = tech.var.prop[set.members], 
                      gene.weights = set.weights, center = center, 
                      gamma = gamma)
    
    results$distance.sq[, i] = vam.results$distance.sq
    results$cdf.value[, i] = vam.results$cdf.value
    pb$tick()
  }
  return(results)
}
getTechVarPropV5 <- function (seurat.data) 
{
  if (!is(seurat.data@assays$RNA, "Assay5")) {
    stop("Trying to call Assay5 logic on non-V5 object!")
  }
  p = nrow(seurat.data@assays$RNA@layers$data)
  meta.data.cols = colnames(seurat.data@assays$RNA@meta.data)
  if (length(meta.data.cols) <= 3 | meta.data.cols[2] != "vf_vst_counts_variance" | 
      meta.data.cols[3] != "vf_vst_counts_variance.expected") {
    message("Did not find vst variance decomposition, setting technical variance proportion to 1")
    return(rep(1, p))
  }
  tech.var.prop = seurat.data@assays$RNA@meta.data[, 3]/seurat.data@assays$RNA@meta.data[, 
                                                                                         2]
  tech.var.prop[which(is.nan(tech.var.prop))] = 1
  return(tech.var.prop)
}

gene_nums <-c(30,50,75,100,150,200,250,300,400,500)

for (gene_num in gene_nums) {
  MOR_genes <- read.csv("HCPMMP/MOR/PLS_report_MOR.csv",row.names=1)
  MOR_genes <- MOR_genes[MOR_genes$Sign>0,]
  MOR_genes <- MOR_genes[1:gene_num,]
  KOR_genes <- read.csv("HCPMMP/KOR/PLS_report_KOR.csv",row.names=1)
  KOR_genes <- KOR_genes[KOR_genes$Sign>0,]
  KOR_genes <- KOR_genes[1:gene_num,]
  
  signatures <- list(
    MOR = rownames(MOR_genes),
    KOR = rownames(KOR_genes)
  )
  # UCell
  seurat_obj <- AddModuleScore_UCell(seurat_obj, features=signatures, name=NULL)
  UCell_scores <- seurat_obj@meta.data[, c("MOR", "KOR"), drop = FALSE]
  write.csv(UCell_scores,file=paste0('Benchmark/Opoid_system/UCell_',gene_num,'genes.csv'))
  # VAM
  VAM_scores <- VAM_Seurat(seurat.data=seurat_obj,
                           gene.set.collection=signatures)
  VAM_scores <- as.data.frame(VAM_scores[["distance.sq"]])
  write.csv(VAM_scores,file=paste0('Benchmark/Opoid_system/VAM_',gene_num,'genes.csv'))
}